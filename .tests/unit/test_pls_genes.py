"""Tests of the step 8 PLS gene assignment."""

import pls_genes

CCRES = [
    ("chr1", 1000, 1200, "EH38D1", "EH38E1", "PLS"),
    ("chr1", 5000, 5150, "EH38D2", "EH38E2", "PLS"),
    ("chr1", 9000, 9200, "EH38D3", "EH38E3", "pELS"),
]
# 0-based base, strand and gene of each TSS on chr1.
TSS = [
    (874, "-", "far_left"),
    (1000, "+", "first_base"),
    (1199, "-", "last_base"),
    (1300, "+", "center_plus_200"),
    (1301, "+", "center_plus_201"),
    (4875, "+", "center_minus_200"),
    (4874, "+", "center_minus_201"),
    (5250, "-", "center_plus_175"),
    (5300, "+", "center_plus_225"),
    (9100, "+", "on_pels"),
]


def write(path, rows):
    """Write rows as tab-separated lines.

    Parameters
    ----------
    path : pathlib.Path
        File to write.
    rows : list of tuple
        Fields per line.
    """
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))


def test_every_tss_within_200_bp_of_the_center(tmp_path):
    """Give a PLS every TSS inside it or within 200 bp of its center."""
    write(tmp_path / "ccres.bed", CCRES)
    write(
        tmp_path / "tss.bed",
        [("chr1", p + 1, p + 1, f"T{i}", ".", s, g) for i, (p, s, g) in enumerate(TSS)],
    )
    pls_genes.main(
        [
            "--ccres",
            str(tmp_path / "ccres.bed"),
            "--tss",
            str(tmp_path / "tss.bed"),
            "--output",
            str(tmp_path / "genes.txt"),
        ]
    )
    assert (tmp_path / "genes.txt").read_text().splitlines() == [
        "EH38E1\tcenter_plus_200",
        "EH38E1\tfirst_base",
        "EH38E1\tlast_base",
        "EH38E2\tcenter_minus_200",
        "EH38E2\tcenter_plus_175",
    ]
