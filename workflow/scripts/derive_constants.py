"""Derive the step 2 DHS signal cutoff from the data.

The release hardcodes dhs_signal_cutoff, and the run keeps that value. This
reports what the data give next to it: the 10th percentile of the DHS-All
signal. The release took 10% of 10,000,000 randomly selected DHSs and kept no
seed.
"""

import argparse

import polars as pl


def build_parser():
    """Build the derive_constants.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the derive_constants.py flags.
    """
    p = argparse.ArgumentParser(description="Derive the DHS signal cutoff.")
    p.add_argument("--dhs-all", required=True, help="DHS-All Parquet of step 2.")
    p.add_argument("--cutoff", required=True, help="Release dhs_signal_cutoff.")
    p.add_argument("-o", "--output", required=True, help="TSV to write.")
    return p


def main(argv=None):
    """Write the cutoff's release value next to its derived value.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    signal = pl.read_parquet(args.dhs_all, columns=["signal"]).to_series()
    cutoff = signal.quantile(0.1, interpolation="linear")
    with open(args.output, "w") as f:
        f.write("constant\trelease\tderived\tderivation\n")
        f.write(
            f"dhs_signal_cutoff\t{args.cutoff}\t{cutoff:.6g}\t"
            "10th percentile of the DHS-All signal\n"
        )


if __name__ == "__main__":
    main()
