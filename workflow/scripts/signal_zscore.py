"""Step 4b: log z-scores of one experiment's signal over the widened anchors.

Port of Toolkit/Retrieve-Signal.sh with log-zscore-normalization.py.
bigWigAverageOverBed averages the signal over each anchor. Anchors with zero
mean0 or zero coverage get z = -inf, the log of 0, and signal 0. The original
gave them z = -10, above any nonzero signal with a lower z. The others get
z = (ln(mean0) - mean) / sd over the nonzero anchors. The original took log10,
which gives the same z up to rounding. Rows are ranked by z, highest first, and
written in anchor name order: name, z, signal and rank.

The mean and SD sum exactly with `math.fsum`, so z does not depend on row
order. The log is glibc's, through polars.
"""

import argparse
import math

import _bed
import polars as pl


def build_parser():
    """Build the signal_zscore.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the signal_zscore.py flags.
    """
    p = argparse.ArgumentParser(description="Log z-score one signal track.")
    p.add_argument("--anchors", required=True, help="Widened anchor Parquet.")
    p.add_argument("--bigwig", required=True, help="Signal bigWig.")
    p.add_argument("-o", "--output", required=True, help="Z-score Parquet to write.")
    return p


def score(out):
    """Score the anchors of one bigWigAverageOverBed table.

    Parameters
    ----------
    out : polars.DataFrame
        The tool's columns: name, size, covered, sum, mean0 and mean.

    Returns
    -------
    polars.DataFrame
        Name, z, signal and rank of each anchor, in name order. Rank is 1 +
        the number of larger z-scores.
    """
    zero = (pl.col("mean0") == 0) | (pl.col("covered") == 0)
    logs = out.filter(~zero)["mean0"].log()
    mean = math.fsum(logs.to_list()) / len(logs)
    deviation = logs - mean
    sd = math.sqrt(math.fsum((deviation * deviation).to_list()) / len(logs))
    if not (deviation / sd).is_finite().all():
        raise ValueError("Non-finite z-score.")
    table = out.select(
        "name",
        z=pl.when(zero)
        .then(float("-inf"))
        .otherwise((pl.col("mean0").log() - mean) / sd),
        signal=pl.when(zero).then(0.0).otherwise("mean0"),
    ).with_columns(rank=pl.col("z").rank("min", descending=True).cast(pl.Int64))
    return table.sort("name")


def main(argv=None):
    """Write name, z, signal and rank of every widened anchor.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    anchors = pl.read_parquet(args.anchors).rename({"accession": "name"})
    _bed.write_parquet(score(_bed.bigwig_averages(anchors, args.bigwig)), args.output)


if __name__ == "__main__":
    main()
