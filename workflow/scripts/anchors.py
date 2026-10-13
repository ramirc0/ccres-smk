"""Anchors: rDHSs and rTFBSs together, sorted with `sort -k1,1 -k2,2n`.

The original repository has no script for this step. Its output equals
`cat rDHS-Filtered.bed rTFBS.bed | sort -k1,1 -k2,2n`.
"""

import argparse

import _bed
import polars as pl

COLUMNS = ["chrom", "start", "end", "accession"]


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
    p.add_argument("-o", "--output", required=True, help="Anchor BED to write.")
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
        [_bed.read_tsv(args.rdhs, COLUMNS), _bed.read_tsv(args.rtfbs, COLUMNS)]
    )
    order = _bed.gnu_sort(anchors, text_keys=["chrom"], numeric_keys=["start"])
    _bed.write_tsv(anchors[order], args.output)


if __name__ == "__main__":
    main()
