"""Tests of the step 5 max z-scores."""

import max_zscore


def test_na_rows_are_skipped(tmp_path):
    """Read no signal file for a mark list row whose file is NA."""
    (tmp_path / "list.txt").write_text(
        "ENCSR000AAA\tENCFF000AAA\tcell_a\n"
        "ENCSR000AAB\tNA\tcell_b\n"
        "ENCSR000AAC\tENCFF000AAC\tcell_c\n"
    )
    signal = tmp_path / "signal"
    signal.mkdir()
    (signal / "ENCSR000AAA-ENCFF000AAA.txt").write_text(
        "EH38D1\t1.5\t2.0\t1\nEH38D2\t-0.5\t1.0\t2\n"
    )
    (signal / "ENCSR000AAC-ENCFF000AAC.txt").write_text(
        "EH38D1\t0.25\t1.0\t2\nEH38D2\t3.0\t9.0\t1\n"
    )
    max_zscore.main(
        [
            "--mark-list",
            str(tmp_path / "list.txt"),
            "--signal",
            str(signal),
            "--output",
            str(tmp_path / "maxZ.txt"),
        ]
    )
    assert (tmp_path / "maxZ.txt").read_text() == "EH38D1\t1.5\nEH38D2\t3.0\n"
