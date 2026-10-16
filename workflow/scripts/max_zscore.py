"""Step 5: the highest z-score of each anchor across one mark's experiments.

Port of 5_Determine-Max-Zscores.sh with select-max-zscore.py. Rows whose file
is NA have no signal file. Every file lists the same anchors in name
order. An anchor's value is its largest z. The original printed it with str().
"""

import argparse
from pathlib import Path

import _bed
import polars as pl


def build_parser():
    """Build the max_zscore.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the max_zscore.py flags.
    """
    p = argparse.ArgumentParser(description="Max z-score per anchor for one mark.")
    p.add_argument("--mark-list", required=True, help="Experiment list of the mark.")
    p.add_argument("--signal", required=True, help="Directory of <exp>-<file>.parquet.")
    p.add_argument("-o", "--output", required=True, help="maxZ Parquet to write.")
    return p


def main(argv=None):
    """Write each anchor's maximum z-score.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.mark_list) as f:
        rows = [line.split("\t") for line in f if line.strip()]
    paths = [
        Path(args.signal, f"{exp}-{acc}.parquet")
        for exp, acc, *_ in rows
        if acc != "NA"
    ]
    first = pl.read_parquet(paths[0], columns=["name", "z"])
    best = first["z"]
    for path in paths[1:]:
        table = pl.read_parquet(path, columns=["name", "z"])
        if not table["name"].equals(first["name"]):
            raise ValueError(f"{path} lists different anchors.")
        best = pl.select(pl.max_horizontal(best, table["z"])).to_series()
    _bed.write_parquet(first.with_columns(z=best), args.output)


if __name__ == "__main__":
    main()
