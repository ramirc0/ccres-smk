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


def test_a_mark_with_no_high_line_still_splits_the_rest():
    """Classify by the next mark when no line of a present mark is high."""
    groups = {"tss": ["D1"], "prox": [], "distal": []}
    classes = cell_type_ccres.dnase_classes(
        groups, [("D1", "0.5")], [("D1", "2.0")], [], [("D1", "3.0")]
    )
    assert classes["pELS"] == ["D1"]
    assert classes["CA-only"] == []


def test_no_dnase_keeps_cres_without_high_marks():
    """Write unclassified rows when no present mark is high in any cCRE."""
    master = [
        ["chr1", "100", "300", "D1", "E1", "PLS"],
        ["chr1", "500", "700", "D2", "E2", "dELS"],
    ]
    groups = {"tss": ["D1"], "prox": [], "distal": ["D2"]}
    lines = cell_type_ccres.no_dnase_lines(
        master, groups, [("D1", "0.5"), ("D2", "0.1")], [], [], [False, True, True]
    )
    assert lines == [
        f"chr1\t{s}\t{e}\t{c}\t0\t.\t{s}\t{e}\t140,140,140\tUnclassified"
        for s, e, c in (("100", "300", "E1"), ("500", "700", "E2"))
    ]
