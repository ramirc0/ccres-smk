"""Tests of the step 4 log z-scores."""

import math

import numpy as np
import polars as pl
import signal_zscore


def averages(mean0, covered):
    """Build a bigWigAverageOverBed table.

    Parameters
    ----------
    mean0 : list of float
        mean0 of each anchor.
    covered : list of int
        Covered bases of each anchor.

    Returns
    -------
    polars.DataFrame
        The tool's columns, anchors named EH38D0000000 onward.
    """
    n = len(mean0)
    return pl.DataFrame(
        {
            "name": [f"EH38D{i:07d}" for i in range(n)],
            "size": [100] * n,
            "covered": covered,
            "sum": mean0,
            "mean0": mean0,
            "mean": mean0,
        }
    )


def test_zero_signal_ranks_last():
    """Give zero-signal anchors z = -inf, below any nonzero signal."""
    table = signal_zscore.score(
        averages([1.0] * 100 + [2.0] * 100 + [1e-30, 0.0], [100] * 201 + [0])
    )
    rows = dict(zip(table["name"], zip(table["z"], table["rank"], table["signal"])))
    assert {rows[f"EH38D{i:07d}"][1] for i in range(100, 200)} == {1}
    assert {rows[f"EH38D{i:07d}"][1] for i in range(100)} == {101}
    assert rows["EH38D0000200"][1] == 201
    assert rows["EH38D0000201"] == (float("-inf"), 202, 0.0)


def test_z_ignores_row_order():
    """Sum exactly, so a permutation of the tool's rows leaves z unchanged."""
    rng = np.random.default_rng(0)
    out = averages(rng.lognormal(0, 3, 50_000).tolist(), [100] * 50_000)
    table = signal_zscore.score(out)
    shuffled = out.sample(fraction=1, shuffle=True, seed=1)
    assert table.equals(signal_zscore.score(shuffled))


def test_log_is_glibc_ln():
    """Take the natural log as Python's math.log does, through glibc."""
    rng = np.random.default_rng(2)
    x = rng.lognormal(0, 5, 200_000)
    expected = np.array([math.log(v) for v in x.tolist()])
    assert np.array_equal(pl.Series(x).log().to_numpy(), expected)
