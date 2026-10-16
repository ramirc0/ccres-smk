"""Step 1: average one experiment's DNase signal over each of its DHSs.

Port of 1_Process-DHSs.sh. DHSs are named `<peak>-<row>` and averaged with
bigWigAverageOverBed. Its mean0 joins them by name. The original sorted both
tables by name and pasted them by position.
"""

import argparse

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
    p.add_argument("--dhs", required=True, help="DHS Parquet from call_dhs.py.")
    p.add_argument(
        "--peak", required=True, help="Peak file accession, the name prefix."
    )
    p.add_argument("--bigwig", required=True, help="DNase signal bigWig.")
    p.add_argument("-o", "--output", required=True, help="Parquet to write.")
    return p


def main(argv=None):
    """Write chrom, start, end, name, FDR and mean0 of each DHS, in calling order.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    dhs = pl.read_parquet(args.dhs).with_columns(
        name=pl.format("{}-{}", pl.lit(args.peak), pl.int_range(1, pl.len() + 1))
    )
    signal = _bed.bigwig_averages(dhs, args.bigwig).select("name", signal="mean0")
    if signal.height != dhs.height:
        raise ValueError(f"{signal.height} averages for {dhs.height} DHSs.")
    processed = dhs.join(signal, on="name", maintain_order="left", validate="1:1")
    _bed.write_parquet(
        processed.select(*_bed.COORDS, "name", "fdr", "signal"), args.output
    )


if __name__ == "__main__":
    main()
