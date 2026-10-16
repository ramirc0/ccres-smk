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
import subprocess
import tempfile
from pathlib import Path

import _bed
import numpy as np
import polars as pl


def build_parser():
    """Build the signal_zscore.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the signal_zscore.py flags.
    """
    p = argparse.ArgumentParser(description="Log z-score one signal track.")
    p.add_argument("--anchors", required=True, help="Widened anchor BED.")
    p.add_argument("--bigwig", required=True, help="Signal bigWig.")
    p.add_argument("-o", "--output", required=True, help="Z-score table to write.")
    return p


def ranks(z):
    """Rank values, highest first. Tied values share a rank.

    Parameters
    ----------
    z : numpy.ndarray
        Values as float64.

    Returns
    -------
    numpy.ndarray
        Rank of each value, 1 + the number of larger values, int64.
    """
    return np.searchsorted(np.sort(-z), -z, side="left") + 1


def score(out):
    """Score the anchors of one bigWigAverageOverBed table.

    Parameters
    ----------
    out : polars.DataFrame
        The tool's string columns: name, size, covered, sum, mean0 and mean.

    Returns
    -------
    polars.DataFrame
        Name, z, signal and rank of each anchor, in anchor name order.
    """
    mean0 = _bed.python_floats(out["mean0"])
    zero = (mean0 == 0) | (_bed.python_floats(out["covered"]) == 0)
    logs = pl.Series(mean0[~zero]).log().to_numpy()
    mean = math.fsum(logs.tolist()) / len(logs)
    deviation = logs - mean
    sd = math.sqrt(math.fsum((deviation * deviation).tolist()) / len(logs))
    z = np.full(len(mean0), -np.inf)
    z[~zero] = (logs - mean) / sd
    if not np.all(np.isfinite(z[~zero])):
        raise ValueError("Non-finite z-score.")

    z_text = [repr(v) for v in z.tolist()]
    signal_text = [
        "0" if is_zero else repr(v) for v, is_zero in zip(mean0.tolist(), zero)
    ]
    table = pl.DataFrame(
        {
            "name": out["name"],
            "z": z_text,
            "signal": signal_text,
            "rank": ranks(z).astype(str).tolist(),
        }
    )
    return table[_bed.gnu_sort(table, ["name"])]


def main(argv=None):
    """Write name, z, signal and rank of every widened anchor.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with tempfile.TemporaryDirectory() as tmp:
        tab = Path(tmp, "out2")
        subprocess.run(
            ["bigWigAverageOverBed", args.bigwig, args.anchors, tab], check=True
        )
        out = _bed.read_tsv(tab, ["name", "size", "covered", "sum", "mean0", "mean"])
    _bed.write_tsv(score(out), args.output)


if __name__ == "__main__":
    main()
