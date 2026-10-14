"""Step 6: classify anchors into cell type-agnostic cCREs and accession them.

Port of 6_Classify-cCREs.sh. Anchors with DNase max z > 1.64 split into TSS,
proximal and distal groups (see `_cres.tss_groups`). Each group then splits by
max z > 1.64 of H3K4me3, H3K27ac and CTCF into PLS, pELS, dELS, CA-H3K4me3,
CA-CTCF and CA. CA anchors overlapping a TF cluster become CA-TF. Anchors with
DNase max z <= 1.64 that overlap a TF cluster become TF, whatever their
H3K4me3, H3K27ac and CTCF z-scores. The paper describes TF cCREs as low in
histone marks too. The cCREs are sorted with `sort -k1,1 -k2,2n` and
accessioned, reusing the accessions of identical V1-V3 cCREs.
"""

import argparse

import _bed
import _cres
import _registry
import polars as pl
from _cres import above, at_most, awk_in

ANCHOR = ["chrom", "start", "end", "accession"]


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


def names(lines):
    """Return the names of signal lines.

    Parameters
    ----------
    lines : list of tuple
        Signal lines, name first.

    Returns
    -------
    list of str
        Names in order.
    """
    return [line[0] for line in lines]


def classify(groups, k4, k27, ctcf):
    """Split the TSS, proximal and distal anchors into classes.

    Parameters
    ----------
    groups : dict of str to list of str
        Anchor names of the tss, prox and distal groups.
    k4, k27, ctcf : list of tuple of (str, str)
        H3K4me3, H3K27ac and CTCF maxZ lines.

    Returns
    -------
    dict of str to list of tuple
        The maxZ lines per class: PLS, pELS, dELS, CA-H3K4me3, CA-CTCF, CA.
    """
    tss_k4 = awk_in(groups["tss"], k4)
    no1 = at_most(tss_k4)
    no2 = at_most(awk_in(names(no1), k27))
    classes = {
        "PLS": above(tss_k4),
        "pELS": above(awk_in(names(no1), k27)),
        "dELS": [],
        "CA-H3K4me3": [],
        "CA-CTCF": above(awk_in(names(no2), ctcf)),
        "CA": at_most(awk_in(names(no2), ctcf)),
    }
    for group, enhancer in (("prox", "pELS"), ("distal", "dELS")):
        group_k27 = awk_in(groups[group], k27)
        classes[enhancer] += above(group_k27)
        no1 = at_most(group_k27)
        classes["CA-H3K4me3"] += above(awk_in(names(no1), k4))
        no2 = at_most(awk_in(names(no1), k4))
        classes["CA-CTCF"] += above(awk_in(names(no2), ctcf))
        classes["CA"] += at_most(awk_in(names(no2), ctcf))
    return classes


def main(argv=None):
    """Write the accessioned cCREs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    anchors = _bed.read_tsv(args.anchors, ANCHOR)
    tss = _bed.read_tsv(args.tss, _cres.TSS)
    window = _bed.read_tsv(args.tss_4k, _cres.TSS)
    tf = _bed.read_tsv(args.tf_clusters, [f"c{i}" for i in range(1, 9)])
    codes = _bed.shared_codes(anchors["chrom"], tss["chrom"], window["chrom"], tf["c1"])
    anchor_iv, tss_iv, window_iv, tf_iv = (
        (c, t[s].cast(pl.Int64).to_numpy(), t[e].cast(pl.Int64).to_numpy())
        for c, t, s, e in zip(
            codes,
            (anchors, tss, window, tf),
            ("start", "start", "start", "c2"),
            ("end", "end", "end", "c3"),
        )
    )
    rows = list(zip(anchors["accession"].to_list(), range(anchors.height)))

    bed = [i for _, i in awk_in(names(above(_cres.read_signal(args.dnase))), rows)]
    groups = _cres.tss_groups(
        [rows[i][0] for i in bed],
        tss,
        tss_iv,
        window_iv,
        tuple(a[bed] for a in anchor_iv),
    )
    classes = classify(
        groups,
        _cres.read_signal(args.h3k4me3),
        _cres.read_signal(args.h3k27ac),
        _cres.read_signal(args.ctcf),
    )

    labelled = []
    for cls in ("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF"):
        labelled += [(i, cls) for _, i in awk_in(names(classes[cls]), rows)]
    on_tf = _bed.overlaps_any(anchor_iv, tf_iv)
    ca = [i for _, i in awk_in(names(classes["CA"]), rows)]
    labelled += [(i, "CA-TF") for i in ca if on_tf[i]]
    labelled += [(i, "CA") for i in ca if not on_tf[i]]
    tf_rows = [rows[i] for i in range(anchors.height) if on_tf[i]]
    labelled_names = [rows[i][0] for i, _ in labelled]
    labelled += [(i, "TF") for _, i in awk_in(labelled_names, tf_rows, negate=True)]

    index, label = zip(*labelled)
    ccres = anchors[list(index)].with_columns(label=pl.Series(label, dtype=pl.String))
    ccres = ccres[_bed.gnu_sort(ccres, text_keys=["chrom"], numeric_keys=["start"])]
    previous = _bed.read_tsv(
        args.previous, ["chrom", "start", "end", "acc", "class", "version"]
    )
    known = _registry.previous_accessions(previous.iter_rows(), "cCRE")
    regions = ccres.select("chrom", "start", "end").iter_rows()
    accession = _registry.accessions(regions, args.genome, "cCRE", known)
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


if __name__ == "__main__":
    main()
