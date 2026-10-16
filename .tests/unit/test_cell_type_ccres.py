"""Tests of the step 7 cell type-specific classification."""

import cell_type_ccres
import polars as pl


def ccres(group, **z):
    """Build one cCRE per group with its z-scores.

    Parameters
    ----------
    group : list of str
        Group of each cCRE.
    **z : list of float or None
        z-scores per mark column. A mark left out has no data.

    Returns
    -------
    polars.DataFrame
        cCREs with chrom, start, end, group and one column per mark.
    """
    n = len(group)
    marks = {m: z.get(m, [None] * n) for m in ("dnase", "h3k4me3", "h3k27ac", "ctcf")}
    return pl.DataFrame(
        {
            "chrom": ["chr1"] * n,
            "start": [1000 * i for i in range(n)],
            "end": [1000 * i + 200 for i in range(n)],
            "group": group,
            **marks,
        },
        schema_overrides=dict.fromkeys(marks, pl.Float64),
    )


def test_tf_summit_is_the_offset_base(tmp_path):
    """Place a summit at start + column 10, the base step 3 uses."""
    path = tmp_path / "peaks.bed"
    path.write_text("chr1\t1000\t1200\t.\t0\t.\t5\t-1\t-1\t50\n")
    assert cell_type_ccres.tf_summits([path]).rows() == [("chr1", 1050, 1051)]


def test_a_mark_with_no_high_line_still_splits_the_rest():
    """Classify by the next mark when no row of a present mark is high."""
    table = ccres(
        ["tss", "distal"], dnase=[2.0, 0.5], h3k4me3=[0.5, 0.5], ctcf=[3.0, 3.0]
    )
    assert cell_type_ccres.dnase_labels(table, []).to_list() == ["CA-CTCF", "Low-DNase"]


def test_no_dnase_keeps_cres_without_high_marks():
    """Write unclassified rows when no present mark is high in any cCRE."""
    colors = cell_type_ccres.no_dnase_colors(
        ccres(["tss", "distal"], h3k4me3=[0.5, 0.1])
    )
    assert colors.rows() == [("140,140,140", "Unclassified")] * 2


def test_no_dnase_label_lists_every_high_mark():
    """Color by the first high mark and list all high marks in the label."""
    colors = cell_type_ccres.no_dnase_colors(
        ccres(["tss", "prox"], h3k4me3=[2.0, 2.0], h3k27ac=[2.0, 0.5], ctcf=[2.0, 2.0])
    )
    assert colors.rows() == [
        ("255,0,0", "High-H3K4me3,High-H3K27ac,High-CTCF"),
        ("255,170,170", "High-H3K4me3,High-CTCF"),
    ]
