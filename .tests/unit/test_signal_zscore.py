"""Tests of the step 4 log z-scores."""

import math

import numpy as np
import polars as pl
import signal_zscore


def test_zero_signal_ranks_last():
    """Give zero-signal anchors z = -inf, below any nonzero signal."""
    names = [f"EH38D{i:07d}" for i in range(202)]
    mean0 = ["1"] * 200 + ["1e-30", "0"]
    out = pl.DataFrame(
        {
            "name": names,
            "size": ["100"] * 202,
            "covered": ["100"] * 201 + ["0"],
            "sum": mean0,
            "mean0": mean0,
            "mean": mean0,
        }
    )
    table = signal_zscore.score(out)
    rows = dict(zip(table["name"], zip(table["z"], table["rank"])))
    assert {rows[n][1] for n in names[:200]} == {"1"}
    assert rows[names[200]][1] == "201"
    assert rows[names[201]] == ("-inf", "202")


def test_z_ignores_row_order():
    """Sum exactly, so a permutation of the tool's rows leaves z unchanged."""
    rng = np.random.default_rng(0)
    mean0 = [repr(v) for v in rng.lognormal(0, 3, 50_000).tolist()]
    out = pl.DataFrame({"name": [f"EH38D{i:07d}" for i in range(50_000)]})
    out = out.with_columns(
        size=pl.lit("100"), covered=pl.lit("100"), sum=pl.Series(mean0)
    ).with_columns(mean0=pl.col("sum"), mean=pl.col("sum"))
    table = signal_zscore.score(out)
    assert table.equals(
        signal_zscore.score(out.sample(fraction=1, shuffle=True, seed=1))
    )


def test_log_is_glibc_ln():
    """Take the natural log as Python's math.log does, through glibc."""
    rng = np.random.default_rng(2)
    x = rng.lognormal(0, 5, 200_000)
    expected = np.array([math.log(v) for v in x.tolist()])
    assert np.array_equal(pl.Series(x).log().to_numpy(), expected)
