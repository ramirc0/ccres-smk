"""Step 4b: log z-scores of one experiment's signal over the widened anchors.

Port of Toolkit/Retrieve-Signal.sh with log-zscore-normalization.py.
bigWigAverageOverBed averages the signal over each anchor. Anchors with zero
mean0 or zero coverage get z = -10 and signal 0. The others get
z = (log10(mean0) - mean) / sd over the nonzero anchors, summed in the tool's
output order. Rows are ranked by z, highest first, and written in anchor
name order: name, z, signal and rank.

`_zscore` reproduces the glibc log and numpy sums of the original run, so the
z-scores match it bit for bit on any host.
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import _bed
import _zscore
import numpy as np
import polars as pl

ZERO_Z = "-10"


def build_parser():
    """Build the signal_zscore.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the signal_zscore.py flags.
    """
    p = argparse.ArgumentParser(description="Log z-score one signal track.")
    p.add_argument("--anchors", required=True, help="Widened anchor BED.")
    p.add_argument(
        "--bigwig", help="Signal bigWig. Without it the output is empty, as for NA."
    )
    p.add_argument("-o", "--output", required=True, help="Z-score table to write.")
    return p


def awk_ranks(z):
    """Rank values like the original's `sort -k2,2rg | awk` rank loop.

    The loop starts with before=0 and takes a new rank only when the value
    differs from the previous one. Ties share the rank of their first row, so
    a value's rank is 1 + the number of larger values. A top value of exactly
    0 would keep rank 0.

    Parameters
    ----------
    z : numpy.ndarray
        Values as float64.

    Returns
    -------
    numpy.ndarray
        Rank of each value, int64.
    """
    descending = -np.sort(-z)
    rank = np.searchsorted(-descending, -z, side="left") + 1
    rank[z == descending[0]] = 0 if descending[0] == 0 else 1
    return rank


def main(argv=None):
    """Write name, z, signal and rank of every widened anchor.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    if args.bigwig is None:
        Path(args.output).write_text("")
        return

    with tempfile.TemporaryDirectory() as tmp:
        tab = Path(tmp, "out2")
        subprocess.run(
            ["bigWigAverageOverBed", args.bigwig, args.anchors, tab], check=True
        )
        out = _bed.read_tsv(tab, ["name", "size", "covered", "sum", "mean0", "mean"])

    mean0 = _bed.python_floats(out["mean0"])
    zero = (mean0 == 0) | (_bed.python_floats(out["covered"]) == 0)
    logs = _zscore.log10(mean0[~zero])
    mean, sd = _zscore.mean_std(logs)
    z = np.full(len(mean0), -10.0)
    z[~zero] = (logs - mean) / sd
    if not np.all(np.isfinite(z)):
        raise ValueError("Non-finite z-score.")

    z_text = [ZERO_Z if is_zero else repr(v) for v, is_zero in zip(z.tolist(), zero)]
    signal_text = [
        "0" if is_zero else repr(v) for v, is_zero in zip(mean0.tolist(), zero)
    ]
    table = pl.DataFrame(
        {
            "name": out["name"],
            "z": z_text,
            "signal": signal_text,
            "rank": awk_ranks(z).astype(str).tolist(),
        }
    )
    _bed.write_tsv(table[_bed.gnu_sort(table, ["name"])], args.output)
    _bed.sync(args.output)


if __name__ == "__main__":
    main()
