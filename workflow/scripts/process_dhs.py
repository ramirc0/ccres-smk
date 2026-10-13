"""Step 1: average one experiment's DNase signal over each of its DHSs.

Port of 1_Process-DHSs.sh. DHSs are named `<peak>-<line>`, sorted by name with
`sort -k4,4`, and averaged with bigWigAverageOverBed. Its output, sorted with
`sort -k1,1`, is pasted next to them by position. Column 6 is the tool's mean0
text, unchanged.
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import _bed
import polars as pl


def build_parser():
    """Build the process_dhs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the process_dhs.py flags.
    """
    p = argparse.ArgumentParser(description="Average DNase signal over DHSs.")
    p.add_argument("--dhs", required=True, help="DHS BED from call_dhs.py.")
    p.add_argument(
        "--peak", required=True, help="Peak file accession, the name prefix."
    )
    p.add_argument("--bigwig", required=True, help="DNase signal bigWig.")
    p.add_argument("-o", "--output", required=True, help="Processed DHS BED to write.")
    return p


def main(argv=None):
    """Write chrom, start, end, name, FDR and mean0 of each DHS, in name order.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    dhs = _bed.read_tsv(args.dhs, ["chrom", "start", "end", "fdr"])
    named = dhs.select(
        "chrom",
        "start",
        "end",
        pl.format("{}-{}", pl.lit(args.peak), pl.int_range(1, pl.len() + 1)).alias(
            "name"
        ),
        "fdr",
    )
    # A -k4,4 key keeps the tab before field 4.
    named = named[_bed.gnu_sort(named.with_columns(key="\t" + pl.col("name")), ["key"])]
    named = named.select("chrom", "start", "end", "name", "fdr")

    with tempfile.TemporaryDirectory() as tmp:
        bed, tab = Path(tmp, "new.bed"), Path(tmp, "out.tab")
        _bed.write_tsv(named.select("chrom", "start", "end", "name"), bed)
        subprocess.run(["bigWigAverageOverBed", args.bigwig, bed, tab], check=True)
        signal = _bed.read_tsv(tab, ["name", "size", "covered", "sum", "mean0", "mean"])
    signal = signal[_bed.gnu_sort(signal, ["name"])]

    if signal.height != named.height:
        raise ValueError(f"{signal.height} averages for {named.height} DHSs.")
    _bed.write_tsv(named.with_columns(mean0=signal["mean0"]), args.output)


if __name__ == "__main__":
    main()
