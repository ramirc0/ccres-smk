"""Tests of the step 7 cell type-specific classification."""

import cell_type_ccres
import polars as pl


def test_tf_summit_is_the_offset_base():
    """Place a summit at start + column 10, the base step 3 uses."""
    peaks = pl.DataFrame(
        [["chr1", "1000", "1200", ".", "0", ".", "5", "-1", "-1", "50"]],
        schema=cell_type_ccres.NARROWPEAK,
        orient="row",
    )
    summits = cell_type_ccres.summit_intervals(peaks)
    assert summits.rows() == [("chr1", 1050, 1051)]
