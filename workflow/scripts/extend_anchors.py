r"""Step 4a: widen the anchors for one signal width and run them through `sort -u`.

Port of the first pipe of Toolkit/Retrieve-Signal.sh:
`printf "%s\t%.0f\t%.0f\t%s\n", $1, $2 - width, $3 + width, $4`, starts below 0
set to 0, then `sort -u` under en_US.UTF-8. Every experiment of a mark averages
its signal over this one file.
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
    p = argparse.ArgumentParser(description="Widen and sort the anchors.")
    p.add_argument("--anchors", required=True, help="Anchor BED.")
    p.add_argument("--width", type=int, required=True, help="bp added to each side.")
    p.add_argument("-o", "--output", required=True, help="Widened BED to write.")
    return p


def main(argv=None):
    """Write the widened, deduplicated and sorted anchors.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    anchors = _bed.read_tsv(args.anchors, ["chrom", "start", "end", "accession"])
    widened = anchors.with_columns(
        start=(pl.col("start").cast(pl.Int64) - args.width)
        .clip(lower_bound=0)
        .cast(pl.String),
        end=(pl.col("end").cast(pl.Int64) + args.width).cast(pl.String),
    )
    with open(args.output, "w") as f:
        f.writelines(
            f"{line}\n" for line in _bed.sort_unique(_bed.lines(widened).to_list())
        )


if __name__ == "__main__":
    main()
