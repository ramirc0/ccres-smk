"""Tests of the step 4 log z-scores."""

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
