"""Step 7b: classify the cCREs in one biosample.

Port of Toolkit/Split-cCREs.{TF,DNase,No-DNase}.sh. The branch follows the
master cell list row: TF when the biosample has DNase and TF data, DNase when
it has DNase only, No-DNase without DNase. A mark the biosample lacks has no
z-scores, so it is high nowhere. The scripts read it from the placeholder
blank.txt. Every cCRE gets one row, sorted by position and labelled by data
completeness. The result is bgzipped, tabix-indexed and written as a bigBed.

DNase and TF: cCREs with DNase z > 1.64 in the biosample split into TSS,
proximal and distal groups, then into PLS, pELS, dELS, CA-H3K4me3 and CA-CTCF
by the biosample's z-scores. Other DNase-high cCREs are CA-only, or CA-TF if
they overlap a summit of the biosample's TF peaks (TF branch). The rest are
Low-DNase.

No-DNase: all cCREs split into groups and get a color by the first high mark,
H3K4me3 first for TSS cCREs, H3K27ac first for the others, CTCF last. A colored
cCRE's label lists its marks above 1.64. The rest are Unclassified. The
script's empty-file fallbacks could give a cCRE two rows when every cCRE of a
group was high in a mark.
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import _bed
import _cres
import polars as pl
from _cres import high

MASTER = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "rdhs": pl.String,
    "ccre": pl.String,
    "class": pl.String,
}
NARROWPEAK = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    **{f"c{i}": pl.String for i in range(4, 10)},
    "offset": pl.Int64,
}
MARKS = {"h3k4me3": "High-H3K4me3", "h3k27ac": "High-H3K27ac", "ctcf": "High-CTCF"}
# label: color, as the Split scripts print them.
COLORS = {
    "CA-CTCF": "0,176,240",
    "dELS": "255,205,0",
    "pELS": "255,167,0",
    "PLS": "255,0,0",
    "CA-H3K4me3": "255,170,170",
    "CA-TF": "190,40,229",
    "CA-only": "6,218,147",
    "Low-DNase": "225,225,225",
    "Unclassified": "140,140,140",
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
    p.add_argument("--dnase", help="DNase z-scores. Absent: No-DNase branch.")
    for mark in ("H3K4me3", "H3K27ac", "CTCF"):
        p.add_argument(f"--{mark.lower()}", help=f"{mark} z-scores. Absent: none.")
    p.add_argument("--tss", required=True, help="TSS BED, zero-length.")
    p.add_argument("--tss-4k", required=True, help="TSS +-2 kb BED.")
    p.add_argument(
        "--tf-peaks",
        nargs="+",
        default=[],
        help="Peaks of the biosample's TF data. Given: TF branch.",
    )
    p.add_argument("--autosql", required=True, help="autoSql of the bigBed fields.")
    p.add_argument("--chrom-sizes", required=True, help="Chromosome sizes.")
    p.add_argument("-o", "--output", required=True, help="BED.gz to write and index.")
    p.add_argument("--bigbed", required=True, help="bigBed to write.")
    return p


def tf_summits(paths):
    """Read the summit bases of TF peak files.

    The summit is start + column 10, as in step 3. Split-cCREs.TF.sh took the
    base before it.

    Parameters
    ----------
    paths : list of str
        Gzipped narrowPeak files.

    Returns
    -------
    polars.DataFrame
        The chrom, start and end of every summit.
    """
    return pl.concat(
        _bed.read_bed(p, NARROWPEAK).select(
            "chrom",
            start=pl.col("start") + pl.col("offset"),
            end=pl.col("start") + pl.col("offset") + 1,
        )
        for p in paths
    )


def dnase_labels(ccres, tf_peaks):
    """Label cCREs by the DNase and TF branches.

    Parameters
    ----------
    ccres : polars.DataFrame
        The cCREs, with group and z-score columns.
    tf_peaks : list of str
        TF peak files of the biosample, none outside the TF branch.

    Returns
    -------
    polars.Series
        Label of each cCRE: a class, CA-TF, CA-only or Low-DNase.
    """
    on_tf = _bed.overlaps_any(ccres, tf_summits(tf_peaks)) if tf_peaks else False
    ca = pl.when(on_tf).then(pl.lit("CA-TF")).otherwise(pl.lit("CA-only"))
    return ccres.select(
        label=pl.when(high("dnase"))
        .then(_cres.label(ca))
        .otherwise(pl.lit("Low-DNase"))
    ).to_series()


def no_dnase_colors(ccres):
    """Color and label cCREs by the No-DNase branch.

    Parameters
    ----------
    ccres : polars.DataFrame
        The cCREs, with group and z-score columns, null for a mark without data.

    Returns
    -------
    polars.DataFrame
        Color and label of each cCRE.
    """
    tss = pl.col("group") == "tss"
    color = (
        pl.when(tss & high("h3k4me3"))
        .then(pl.lit("255,0,0"))
        .when(high("h3k27ac"))
        .then(pl.lit("255,205,0"))
        .when(~tss & high("h3k4me3"))
        .then(pl.lit("255,170,170"))
        .when(high("ctcf"))
        .then(pl.lit("0,176,240"))
    )
    marks = pl.concat_list(
        [pl.when(high(mark)).then(pl.lit(name)) for mark, name in MARKS.items()]
    ).list.drop_nulls()
    return ccres.select(
        color=color.fill_null(pl.lit(COLORS["Unclassified"])),
        label=pl.when(color.is_null())
        .then(pl.lit("Unclassified"))
        .otherwise(marks.list.join(",")),
    )


def main(argv=None):
    """Write the biosample's classified cCREs as an indexed BED.gz and a bigBed.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    tables = {mark: getattr(args, mark) for mark in ("dnase", *MARKS)}
    ccres = _cres.with_z(_bed.read_bed(args.master, MASTER), "rdhs", tables)
    ccres = ccres.with_columns(
        _cres.tss_groups(
            ccres, _cres.read_tss(args.tss), _bed.read_bed(args.tss_4k, _cres.TSS)
        )
    )
    if args.dnase:
        label = dnase_labels(ccres, args.tf_peaks)
        colored = ccres.with_columns(label, color=label.replace_strict(COLORS))
        complete = all(tables[mark] for mark in MARKS)
    else:
        colored = pl.concat([ccres, no_dnase_colors(ccres)], how="horizontal")
        complete = False

    final = colored.sort(*_bed.COORDS, "ccre").select(
        *_bed.COORDS,
        "ccre",
        score=pl.lit("0"),
        strand=pl.lit("."),
        thick_start="start",
        thick_end="end",
        color="color",
        label="label",
        completeness=pl.lit(COMPLETE if complete else PARTIAL),
    )
    with tempfile.TemporaryDirectory() as tmp:
        bed_path = Path(tmp, "cres.bed")
        _bed.write_bed(final, bed_path)
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


if __name__ == "__main__":
    main()
