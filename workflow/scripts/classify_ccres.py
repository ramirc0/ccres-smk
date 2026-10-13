"""Step 6: classify anchors into cell type-agnostic cCREs and accession them.

Port of 6_Classify-cCREs.sh with calculate-center-distance.py (agnostic).
Anchors with DNase max z > 1.64 split by TSS distance:

- TSS: overlapping a TSS, or proximal with the closest TSS within 200 bp of
  the anchor center.
- proximal: within a TSS +-2 kb window, otherwise.
- distal: outside every window, or proximal with the closest TSS over 2 kb
  from the center.

Each group then splits by max z > 1.64 of H3K4me3, H3K27ac and CTCF into PLS,
pELS, dELS, CA-H3K4me3, CA-CTCF and CA. CA anchors overlapping a TF cluster
become CA-TF. Any other anchor overlapping a TF cluster becomes TF. The cCREs
are sorted with `sort -k1,1 -k2,2n` and accessioned, reusing the accessions of
identical V1-V3 cCREs.
"""

import argparse

import _bed
import _registry
import numpy as np
import polars as pl

HIGH = 1.64
CENTER = 200
FAR = 2000
ANCHOR = ["chrom", "start", "end", "accession"]
TSS = ["chrom", "start", "end", "transcript", "score", "strand", "gene"]


def build_parser():
    """Build the classify_ccres.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the classify_ccres.py flags.
    """
    p = argparse.ArgumentParser(description="Classify cell type-agnostic cCREs.")
    p.add_argument("--anchors", required=True, help="Anchor BED.")
    for mark in ("DNase", "H3K4me3", "H3K27ac", "CTCF"):
        p.add_argument(f"--{mark.lower()}", required=True, help=f"{mark} maxZ table.")
    p.add_argument("--tss", required=True, help="TSS BED, zero-length.")
    p.add_argument("--tss-4k", required=True, help="TSS +-2 kb BED.")
    p.add_argument("--tf-clusters", required=True, help="TF-Clusters.txt.")
    p.add_argument("--previous", required=True, help="V1-V3 cCRE BED.")
    p.add_argument("--genome", required=True, help="hg38 or mm10.")
    p.add_argument("-o", "--output", required=True, help="cCRE BED to write.")
    return p


def high_and_low(path):
    """Split a maxZ table's names by value, as awk's `$2 > 1.64` and `<= 1.64`.

    Parameters
    ----------
    path : str
        The maxZ table of name and value.

    Returns
    -------
    tuple of set of str
        Names above 1.64 and names at or below it. NaN is in neither.
    """
    table = _bed.read_tsv(path, ["name", "value"])
    value = _bed.python_floats(table["value"])
    names = np.array(table["name"].to_list(), dtype=object)
    return set(names[value > HIGH]), set(names[value <= HIGH])


def intervals(df, codes):
    """Return the interval arrays of a table.

    Parameters
    ----------
    df : polars.DataFrame
        Table with start and end columns.
    codes : numpy.ndarray
        Chromosome code of each row.

    Returns
    -------
    tuple of numpy.ndarray
        Codes, starts and ends.
    """
    return (
        codes,
        df["start"].cast(pl.Int64).to_numpy(),
        df["end"].cast(pl.Int64).to_numpy(),
    )


class Regions:
    """Tables placed on one chromosome numbering, for overlap queries.

    Parameters
    ----------
    **tables : polars.DataFrame
        Tables with chrom, start and end columns, by name.
    """

    def __init__(self, **tables):
        self.tables = tables
        codes = _bed.shared_codes(*(t["chrom"] for t in tables.values()))
        self.codes = dict(zip(tables, codes))

    def subset(self, name, rows):
        """Return the interval arrays of some rows of a table, in the given order.

        Parameters
        ----------
        name : str
            Table name.
        rows : numpy.ndarray
            Row indices.

        Returns
        -------
        tuple of numpy.ndarray
            Codes, starts and ends.
        """
        return intervals(self.tables[name][rows], self.codes[name][rows])

    def whole(self, name):
        """Return the interval arrays of a whole table.

        Parameters
        ----------
        name : str
            Table name.

        Returns
        -------
        tuple of numpy.ndarray
            Codes, starts and ends.
        """
        return intervals(self.tables[name], self.codes[name])


def center_distances(prox, tss):
    """Signed distance from each anchor's center to its first closest TSS.

    Port of calculate-center-distance.py (agnostic) on `bedtools closest -d`
    output: the first line per anchor wins, the TSS start minus the center on
    the + strand, the center minus the TSS start on the - strand.

    Parameters
    ----------
    prox : tuple of (polars.DataFrame, tuple of numpy.ndarray)
        Proximal anchors in `sort -k1,1 -k2,2n` order and their intervals.
    tss : tuple of (polars.DataFrame, tuple of numpy.ndarray)
        TSS table and its intervals.

    Returns
    -------
    numpy.ndarray
        Distance per anchor, float64.
    """
    (anchors, a_iv), (sites, b_iv) = prox, tss
    ia, ib = _bed.closest_pairs(a_iv, b_iv)
    first = np.flatnonzero(np.diff(ia, prepend=-1))
    if len(first) != anchors.height:
        raise ValueError("An anchor has no TSS on its chromosome.")
    site = ib[first]
    center = (a_iv[1] + a_iv[2]) / 2
    start = sites["start"].cast(pl.Int64).to_numpy()[site]
    plus = sites["strand"].to_numpy()[site] == "+"
    return np.where(plus, start - center, center - start)


def main(argv=None):
    """Write the accessioned cCREs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    dnase_high, _ = high_and_low(args.dnase)
    k4_high, k4_low = high_and_low(args.h3k4me3)
    k27_high, k27_low = high_and_low(args.h3k27ac)
    ctcf_high, ctcf_low = high_and_low(args.ctcf)

    anchors = _bed.read_tsv(args.anchors, ANCHOR)
    tss = _bed.read_tsv(args.tss, TSS)
    regions = Regions(
        anchors=anchors,
        tss=tss,
        window=_bed.read_tsv(args.tss_4k, TSS),
        tf=_bed.read_tsv(args.tf_clusters, [f"c{i}" for i in range(1, 9)]).rename(
            {"c1": "chrom", "c2": "start", "c3": "end"}
        ),
    )
    name = anchors["accession"].to_numpy()

    def member(names):
        """Flag the anchors whose accession is in `names`.

        Parameters
        ----------
        names : set of str
            Anchor accessions.

        Returns
        -------
        numpy.ndarray
            Flag per anchor.
        """
        return anchors["accession"].is_in(list(names)).to_numpy()

    bed = member(dnase_high)
    on_tss = np.zeros(len(name), dtype=bool)
    on_tss[bed] = _bed.overlaps_any(
        regions.subset("anchors", np.flatnonzero(bed)), regions.whole("tss")
    )
    in_window = _bed.overlaps_any(regions.whole("anchors"), regions.whole("window"))
    proximal = bed & ~on_tss & in_window
    distal = bed & ~in_window

    prox_rows = np.flatnonzero(proximal)
    prox_table = anchors[prox_rows]
    prox_order = _bed.gnu_sort(prox_table, text_keys=["chrom"], numeric_keys=["start"])
    prox_rows = prox_rows[prox_order]
    distance = center_distances(
        (anchors[prox_rows], regions.subset("anchors", prox_rows)),
        (tss, regions.whole("tss")),
    )
    near = (distance >= -CENTER) & (distance <= CENTER)
    far = (distance < -FAR) | (distance > FAR)
    on_tss[prox_rows[near]] = True
    distal[prox_rows[far]] = True
    proximal[prox_rows[near | far]] = False

    groups = {
        "tss": set(name[on_tss]),
        "prox": set(name[proximal]),
        "distal": set(name[distal]),
    }
    classes = classify(groups, k4_high, k4_low, k27_high, k27_low, ctcf_high, ctcf_low)

    on_tf = _bed.overlaps_any(regions.whole("anchors"), regions.whole("tf"))
    label = np.full(len(name), "", dtype=object)
    for cls in ("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF"):
        label[member(classes[cls])] = cls
    ca = member(classes["CA"])
    label[ca & on_tf] = "CA-TF"
    label[ca & ~on_tf] = "CA"
    label[(label == "") & on_tf] = "TF"

    ccres = anchors.with_columns(label=pl.Series(label.tolist(), dtype=pl.String))
    ccres = ccres.filter(pl.col("label") != "")
    ccres = ccres[_bed.gnu_sort(ccres, text_keys=["chrom"], numeric_keys=["start"])]
    previous = _bed.read_tsv(
        args.previous, ["chrom", "start", "end", "acc", "class", "version"]
    )
    known = _registry.previous_accessions(previous.iter_rows(), "cCRE")
    regions_iter = ccres.select("chrom", "start", "end").iter_rows()
    accession = _registry.accessions(regions_iter, args.genome, "cCRE", known)
    _bed.write_tsv(
        ccres.select(
            "chrom",
            "start",
            "end",
            "accession",
            pl.Series("ccre", accession, dtype=pl.String),
            "label",
        ),
        args.output,
    )


def classify(groups, k4_high, k4_low, k27_high, k27_low, ctcf_high, ctcf_low):
    """Split the TSS, proximal and distal anchors into classes.

    Parameters
    ----------
    groups : dict of str to set of str
        Anchor names of the tss, prox and distal groups.
    k4_high, k4_low, k27_high, k27_low, ctcf_high, ctcf_low : set of str
        Names above and at or below 1.64 for H3K4me3, H3K27ac and CTCF.

    Returns
    -------
    dict of str to set of str
        Anchor names per class: PLS, pELS, dELS, CA-H3K4me3, CA-CTCF, CA.
    """
    tss, prox, distal = groups["tss"], groups["prox"], groups["distal"]
    no1 = tss & k4_low
    no2 = no1 & k27_low
    tss_no2 = no2
    prox_no1 = prox & k27_low
    prox_no2 = prox_no1 & k4_low
    distal_no1 = distal & k27_low
    distal_no2 = distal_no1 & k4_low
    rest = tss_no2 | prox_no2 | distal_no2
    return {
        "PLS": tss & k4_high,
        "pELS": (no1 & k27_high) | (prox & k27_high),
        "dELS": distal & k27_high,
        "CA-H3K4me3": (prox_no1 & k4_high) | (distal_no1 & k4_high),
        "CA-CTCF": rest & ctcf_high,
        "CA": rest & ctcf_low,
    }


if __name__ == "__main__":
    main()
