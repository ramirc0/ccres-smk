"""Anchors: rDHSs and rTFBSs together, sorted by position.

The original repository has no script for this step. Its output equals
`cat rDHS-Filtered.bed rTFBS.bed | sort -k1,1 -k2,2n`.
"""

import argparse

import _bed
import polars as pl

ANCHOR = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "accession": pl.String,
}


def build_parser():
    """Build the anchors.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the anchors.py flags.
    """
    p = argparse.ArgumentParser(description="Combine rDHSs and rTFBSs.")
    p.add_argument("--rdhs", required=True, help="rDHS BED.")
    p.add_argument("--rtfbs", required=True, help="rTFBS BED.")
    p.add_argument("-o", "--output", required=True, help="Anchor Parquet to write.")
    return p


def main(argv=None):
    """Write the sorted anchors.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    anchors = pl.concat(
        [_bed.read_bed(args.rdhs, ANCHOR), _bed.read_bed(args.rtfbs, ANCHOR)]
    )
    _bed.write_parquet(anchors.sort(anchors.columns), args.output)


if __name__ == "__main__":
    main()
