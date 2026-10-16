r"""Step 4a: widen the anchors for one signal width and drop repeats.

Port of the first pipe of Toolkit/Retrieve-Signal.sh:
`printf "%s\t%.0f\t%.0f\t%s\n", $1, $2 - width, $3 + width, $4`, starts below 0
set to 0, then `sort -u`. Every experiment of a mark averages its signal over
this one table.
"""

import argparse

import _bed
import polars as pl


def build_parser():
    """Build the extend_anchors.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the extend_anchors.py flags.
    """
    p = argparse.ArgumentParser(description="Widen the anchors.")
    p.add_argument("--anchors", required=True, help="Anchor Parquet.")
    p.add_argument("--width", type=int, required=True, help="bp added to each side.")
    p.add_argument("-o", "--output", required=True, help="Parquet to write.")
    return p


def main(argv=None):
    """Write the widened, deduplicated and sorted anchors.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    widened = pl.read_parquet(args.anchors).with_columns(
        start=(pl.col("start") - args.width).clip(lower_bound=0),
        end=pl.col("end") + args.width,
    )
    _bed.write_parquet(widened.unique().sort(widened.columns), args.output)


if __name__ == "__main__":
    main()
