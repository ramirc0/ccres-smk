"""Tests of the step 5 max z-scores."""

import max_zscore
import polars as pl


def test_na_rows_are_skipped(tmp_path):
    """Read no signal file for a mark list row whose file is NA."""
    (tmp_path / "list.txt").write_text(
        "ENCSR000AAA\tENCFF000AAA\tcell_a\n"
        "ENCSR000AAB\tNA\tcell_b\n"
        "ENCSR000AAC\tENCFF000AAC\tcell_c\n"
    )
    signal = tmp_path / "signal"
    signal.mkdir()
    names = ["EH38D1", "EH38D2"]
    pl.DataFrame({"name": names, "z": [1.5, -0.5]}).write_parquet(
        signal / "ENCSR000AAA-ENCFF000AAA.parquet"
    )
    pl.DataFrame({"name": names, "z": [0.25, 3.0]}).write_parquet(
        signal / "ENCSR000AAC-ENCFF000AAC.parquet"
    )
    out = tmp_path / "maxZ.parquet"
    max_zscore.main(
        [
            "--mark-list",
            str(tmp_path / "list.txt"),
            "--signal",
            str(signal),
            "-o",
            str(out),
        ]
    )
    assert pl.read_parquet(out).rows() == [("EH38D1", 1.5), ("EH38D2", 3.0)]
