"""Step 7b: classify the cCREs in one biosample.

Port of Toolkit/Split-cCREs.{TF,DNase,No-DNase}.sh. The branch follows the
master cell list row: TF when the biosample has DNase and TF data, DNase when
it has DNase only, No-DNase without DNase. A mark the biosample lacks reads as
the placeholder blank.txt. Every cCRE of the release appears once per class it
gets, then rows are sorted with `sort -k1,1 -k2,2n` and labelled by data
completeness. The result is bgzipped, tabix-indexed and written as a bigBed.

DNase and TF: cCREs with DNase z > 1.64 in the biosample split into TSS,
proximal and distal groups, then into PLS, pELS, dELS, CA-H3K4me3 and CA-CTCF
by the biosample's z-scores. Other DNase-high cCREs are CA-only, or CA-TF if
they overlap a summit of the biosample's TF peaks (TF branch). The rest are
Low-DNase.

No-DNase: all cCREs split into groups and get one row per colored class, with
High-<mark> suffixes for the marks above 1.64. The rest are Unclassified.
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import _bed
import _cres
import _encode
import polars as pl
from _cres import above, at_most, awk_in

MASTER = ["chrom", "start", "end", "rdhs", "ccre", "class"]
NARROWPEAK = [f"c{i}" for i in range(1, 11)]
# label: color, as the Split scripts print them.
DNASE_COLORS = {
    "CA-CTCF": "0,176,240",
    "dELS": "255,205,0",
    "pELS": "255,167,0",
    "PLS": "255,0,0",
    "CA-H3K4me3": "255,170,170",
    "CA-TF": "190,40,229",
    "CA-only": "6,218,147",
}
NO_DNASE_COLORS = {
    "blue": "0,176,240",
    "pink": "255,170,170",
    "orange": "255,205,0",
    "red": "255,0,0",
}
COMPLETE = "All-data/Full-classification"
PARTIAL = "Missing-data/Partial-classification"


def build_parser():
    """Build the cell_type_ccres.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the cell_type_ccres.py flags.
    """
    p = argparse.ArgumentParser(description="Classify cCREs in one biosample.")
    p.add_argument("--master", required=True, help="cCREs-Unfiltered BED.")
    p.add_argument("--dnase", help="DNase z-score table. Absent: No-DNase branch.")
    for mark in ("H3K4me3", "H3K27ac", "CTCF"):
        p.add_argument(f"--{mark.lower()}", help=f"{mark} z-scores. Absent: blank.")
    p.add_argument("--tss", required=True, help="TSS BED, zero-length.")
    p.add_argument("--tss-4k", required=True, help="TSS +-2 kb BED.")
    p.add_argument("--tf-list", help="TF experiment list. Given: TF branch.")
    p.add_argument(
        "--biosample", help="Biosample, matched to the TF list's last token."
    )
    p.add_argument("--encode-root", help="Local ENCODE mirror.")
    p.add_argument("--download-dir", help="Downloaded ENCODE files.")
    p.add_argument("--autosql", required=True, help="autoSql of the bigBed fields.")
    p.add_argument("--chrom-sizes", required=True, help="Chromosome sizes.")
    p.add_argument("-o", "--output", required=True, help="BED.gz to write and index.")
    p.add_argument("--bigbed", required=True, help="bigBed to write.")
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


def dnase_classes(groups, k4, k27, ctcf, dnase_high):
    """Split DNase-high cCREs into the classes of the DNase and TF branches.

    Parameters
    ----------
    groups : dict of str to list of str
        The rDHS accessions of the tss, prox and distal groups.
    k4, k27, ctcf : list of tuple of (str, str)
        Signal lines of the biosample, or the blank placeholder line.
    dnase_high : list of tuple of (str, str)
        DNase lines above 1.64.

    Returns
    -------
    dict of str to list of str
        The rDHS accessions per class, CA-only holding the unclassified
        DNase-high ones.
    """
    high = {
        "k4": names(above(k4)),
        "k27": names(above(k27)),
        "ctcf": names(above(ctcf)),
    }
    classes = {"PLS": awk_in(high["k4"], groups["tss"])}
    no1 = awk_in(high["k4"], groups["tss"], negate=True)
    classes["pELS"] = awk_in(high["k27"], no1)
    no2 = awk_in(high["k27"], no1, negate=True)
    classes["CA-CTCF"] = awk_in(high["ctcf"], no2)
    classes["dELS"], classes["CA-H3K4me3"] = [], []
    for group, enhancer in (("prox", "pELS"), ("distal", "dELS")):
        classes[enhancer] += awk_in(high["k27"], groups[group])
        no1 = awk_in(high["k27"], groups[group], negate=True)
        classes["CA-H3K4me3"] += awk_in(high["k4"], no1)
        no2 = awk_in(high["k4"], no1, negate=True)
        classes["CA-CTCF"] += awk_in(high["ctcf"], no2)
    found = [
        n for c in ("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF") for n in classes[c]
    ]
    # An empty all.txt would trip awk_in, so the script wrote "done" into it.
    classes["CA-only"] = awk_in(found or ["done"], names(dnase_high), negate=True)
    return classes


def no_dnase_classes(groups, k4, k27, ctcf):
    """Colored classes of the No-DNase branch, with its fallbacks.

    An empty no1 or no2 falls back to the whole group or to no1, as the script's
    `wc -l` checks do. The fallback can put one cCRE in two classes.

    Parameters
    ----------
    groups : dict of str to list of str
        The rDHS accessions of the tss, prox and distal groups.
    k4, k27, ctcf : list of tuple of (str, str)
        Signal lines of the biosample, or the blank placeholder line.

    Returns
    -------
    dict of str to list of str
        The rDHS accessions per color: red, orange, pink and blue.
    """
    tss_k4 = awk_in(groups["tss"], k4)
    color = {"red": names(above(tss_k4)), "orange": [], "pink": [], "blue": []}
    no1 = names(at_most(tss_k4)) or groups["tss"]
    tss_k27 = awk_in(no1, k27)
    color["orange"] += names(above(tss_k27))
    no2 = names(at_most(tss_k27)) or no1
    color["blue"] += names(above(awk_in(no2, ctcf)))
    for group in ("prox", "distal"):
        group_k27 = awk_in(groups[group], k27)
        color["orange"] += names(above(group_k27))
        no1 = names(at_most(group_k27)) or groups[group]
        no1_k4 = awk_in(no1, k4)
        color["pink"] += names(above(no1_k4))
        no2 = names(at_most(no1_k4)) or no1
        color["blue"] += names(above(awk_in(no2, ctcf)))
    return color


def master_rows(master, accessions):
    """Return master rows whose rDHS is in `accessions`, in master order.

    Parameters
    ----------
    master : list of list of str
        The fields of the cCREs-Unfiltered rows.
    accessions : list of str
        The rDHS accessions, the first file of the awk join.

    Returns
    -------
    list of list of str
        Selected rows.
    """
    return awk_in(accessions, master, field=3)


def class_line(row, color, label):
    """Print a master row as the Split scripts' l.bed line, rDHS first.

    Parameters
    ----------
    row : list of str
        The cCREs-Unfiltered fields.
    color, label : str
        Columns 10 and 11.

    Returns
    -------
    str
        The line: rDHS, chrom, start, end, cCRE, 0, ., start, end, color,
        label.
    """
    chrom, start, end, rdhs, ccre = row[:5]
    fields = (rdhs, chrom, start, end, ccre, "0", ".", start, end, color, label)
    return "\t".join(fields)


def drop_first(line):
    """Drop the first field like `awk '{printf "%s", $2; ...}'`, splitting on blanks.

    Parameters
    ----------
    line : str
        Line of whitespace-separated fields.

    Returns
    -------
    str
        Fields 2 onward, tab-joined. Empty fields vanish, as in awk.
    """
    return "\t".join(line.split()[1:])


def sorted_lines(lines, text_key, numeric_key=None):
    """Sort lines like `sort -k<text>,<text> [-k<num>,<num>n]` under en_US.UTF-8.

    Parameters
    ----------
    lines : list of str
        Tab-separated lines.
    text_key : int
        0-based field compared by collation.
    numeric_key : int, optional
        0-based field compared as an integer.

    Returns
    -------
    list of str
        Sorted lines.
    """
    if not lines:
        return []
    line = pl.Series(lines)
    fields = line.str.split("\t").list
    table = pl.DataFrame(
        {
            "text": fields.get(text_key),
            "line": line,
            **({"number": fields.get(numeric_key)} if numeric_key else {}),
        }
    )
    order = _bed.gnu_sort(
        table,
        text_keys=["text"],
        numeric_keys=["number"] if numeric_key else [],
        line="line",
    )
    return [lines[i] for i in order]


def summit_intervals(peaks):
    """Return the summit base of each narrowPeak row.

    The summit is start + column 10, as in step 3. Split-cCREs.TF.sh took the
    base before it.

    Parameters
    ----------
    peaks : polars.DataFrame
        The narrowPeak string columns c1-c10.

    Returns
    -------
    polars.DataFrame
        The chrom, start and end of each summit.
    """
    summit = pl.col("c2").cast(pl.Int64) + pl.col("c10").cast(pl.Int64)
    return peaks.select(chrom="c1", start=summit, end=summit + 1)


def tf_summits(args):
    """Read the summits of the biosample's TF peaks.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed flags with tf_list, biosample, encode_root and download_dir.

    Returns
    -------
    polars.DataFrame
        The chrom, start and end of every summit.
    """
    with open(args.tf_list) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    parts = []
    for row in rows:
        if row[-1].split()[-1] != args.biosample:
            continue
        path = _encode.encode_path(
            args.encode_root, args.download_dir, row[0], row[1], "bed.gz"
        )
        parts.append(summit_intervals(_bed.read_tsv(path, NARROWPEAK)))
    return (
        pl.concat(parts)
        if parts
        else pl.DataFrame(
            schema={"chrom": pl.String, "start": pl.Int64, "end": pl.Int64}
        )
    )


def main(argv=None):
    """Write the biosample's classified cCREs as an indexed BED.gz and a bigBed.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    master = [list(r) for r in _bed.read_tsv(args.master, MASTER).rows()]
    tss = _bed.read_tsv(args.tss, _cres.TSS)
    window = _bed.read_tsv(args.tss_4k, _cres.TSS)
    k4, k27, ctcf = (
        _cres.read_signal(getattr(args, m)) for m in ("h3k4me3", "h3k27ac", "ctcf")
    )
    blank = [getattr(args, m) is None for m in ("h3k4me3", "h3k27ac", "ctcf")]

    if args.dnase:
        dnase_high = above(_cres.read_signal(args.dnase))
        bed = master_rows(master, names(dnase_high))
    else:
        bed = master
    bed_table = pl.DataFrame(bed, schema=dict.fromkeys(MASTER, pl.String), orient="row")
    codes = _bed.shared_codes(bed_table["chrom"], tss["chrom"], window["chrom"])
    region_iv, tss_iv, window_iv = (
        (c, t["start"].cast(pl.Int64).to_numpy(), t["end"].cast(pl.Int64).to_numpy())
        for c, t in zip(codes, (bed_table, tss, window))
    )
    groups = _cres.tss_groups(
        bed_table["rdhs"].to_list(), tss, tss_iv, window_iv, region_iv
    )

    if args.dnase:
        lines = dnase_lines(args, master, groups, k4, k27, ctcf, dnase_high)
        completeness = PARTIAL if any(blank) else COMPLETE
    else:
        lines = no_dnase_lines(master, groups, k4, k27, ctcf, blank)
        completeness = PARTIAL
    final = sorted_lines(lines, 0, 1)

    with tempfile.TemporaryDirectory() as tmp:
        bed_path = Path(tmp, "cres.bed")
        bed_path.write_text("".join(f"{line}\t{completeness}\n" for line in final))
        with open(args.output, "wb") as out:
            subprocess.run(["bgzip", "-c", bed_path], stdout=out, check=True)
    subprocess.run(["tabix", "-f", "-p", "bed", args.output], check=True)
    subprocess.run(
        [
            "bedToBigBed",
            "-type=bed9+2",
            f"-as={args.autosql}",
            args.output,
            args.chrom_sizes,
            args.bigbed,
        ],
        check=True,
    )


def dnase_lines(args, master, groups, k4, k27, ctcf, dnase_high):
    """Build the l.bed and p.bed lines of the DNase and TF branches.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed flags. A TF list selects the TF branch.
    master : list of list of str
        The cCREs-Unfiltered rows.
    groups : dict of str to list of str
        The rDHS accessions of the tss, prox and distal groups.
    k4, k27, ctcf : list of tuple of (str, str)
        Signal lines of the biosample, or the blank placeholder line.
    dnase_high : list of tuple of (str, str)
        DNase lines above 1.64.

    Returns
    -------
    list of str
        Unsorted lines, rDHS column dropped.
    """
    classes = dnase_classes(groups, k4, k27, ctcf, dnase_high)
    lines = []
    for label in ("CA-CTCF", "dELS", "pELS", "PLS", "CA-H3K4me3"):
        lines += [
            class_line(r, DNASE_COLORS[label], label)
            for r in master_rows(master, classes[label])
        ]
    ca = master_rows(master, classes["CA-only"])
    if args.tf_list:
        summits = tf_summits(args)
        ca_table = pl.DataFrame(
            ca, schema=dict.fromkeys(MASTER, pl.String), orient="row"
        )
        ca_codes, summit_codes = _bed.shared_codes(ca_table["chrom"], summits["chrom"])
        hit = _bed.overlaps_any(
            (
                ca_codes,
                ca_table["start"].cast(pl.Int64).to_numpy(),
                ca_table["end"].cast(pl.Int64).to_numpy(),
            ),
            (summit_codes, summits["start"].to_numpy(), summits["end"].to_numpy()),
        )
        lines += [
            class_line(r, DNASE_COLORS["CA-TF"], "CA-TF") for r, h in zip(ca, hit) if h
        ]
        ca = [r for r, h in zip(ca, hit) if not h]
    lines += [class_line(r, DNASE_COLORS["CA-only"], "CA-only") for r in ca]
    lines = [drop_first(line) for line in lines]
    labelled = [line.split("\t")[3] for line in lines]
    low = awk_in(labelled, master, negate=True, field=4)
    return lines + [
        "\t".join(
            [r[0], r[1], r[2], r[4], "0", ".", r[1], r[2], "225,225,225", "Low-DNase"]
        )
        for r in low
    ]


def no_dnase_lines(master, groups, k4, k27, ctcf, blank):
    """Build the m.bed and p.bed lines of the No-DNase branch.

    Parameters
    ----------
    master : list of list of str
        The cCREs-Unfiltered rows.
    groups : dict of str to list of str
        The rDHS accessions of the tss, prox and distal groups.
    k4, k27, ctcf : list of tuple of (str, str)
        Signal lines of the biosample, or the blank placeholder line.
    blank : list of bool
        Whether H3K4me3, H3K27ac and CTCF are the blank placeholder.

    Returns
    -------
    list of str
        Unsorted lines, rDHS column dropped.
    """
    color = no_dnase_classes(groups, k4, k27, ctcf)
    lines = []
    for name in ("blue", "pink", "orange", "red"):
        lines += [
            class_line(r, NO_DNASE_COLORS[name], "#")
            for r in master_rows(master, color[name])
        ]
    for signal, is_blank, suffix in zip(
        (k4, k27, ctcf), blank, ("High-H3K4me3,", "High-H3K27ac,", "High-CTCF,")
    ):
        if not is_blank:
            first = [line.split("\t")[0] for line in lines]
            high = names(above(awk_in(first, signal)))
            marked = [line + suffix for _, line in awk_in(high, _keyed(lines))]
            rest = [line for _, line in awk_in(high, _keyed(lines), negate=True)]
            lines = marked + rest
        lines = sorted_lines(lines, 0)
    lines = [drop_first(line.removesuffix(",").replace("#", "")) for line in lines]
    labelled = [line.split("\t")[3] for line in lines]
    unclassified = awk_in(labelled, master, negate=True, field=4)
    return lines + [
        "\t".join(
            [
                r[0],
                r[1],
                r[2],
                r[4],
                "0",
                ".",
                r[1],
                r[2],
                "140,140,140",
                "Unclassified",
            ]
        )
        for r in unclassified
    ]


def _keyed(lines):
    """Key tab-separated lines by their first field.

    Parameters
    ----------
    lines : list of str
        Lines.

    Returns
    -------
    list of tuple of (str, str)
        First field and line.
    """
    return [(line.split("\t")[0], line) for line in lines]


if __name__ == "__main__":
    main()
