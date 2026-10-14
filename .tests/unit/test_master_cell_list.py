"""Tests of the step 7 master cell list."""

import master_cell_list as mcl


def write_lists(tmp_path, lists, tf_rows=()):
    """Write the four mark lists and a TF list.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to write in.
    lists : dict of str to list of list of str
        Rows per mark.
    tf_rows : sequence of list of str
        TF list rows.

    Returns
    -------
    tuple of (list of pathlib.Path, pathlib.Path)
        Mark list paths in `MARKS` order, then the TF list path.
    """
    paths = []
    for mark in mcl.MARKS:
        path = tmp_path / f"{mark}-List.txt"
        path.write_text("".join("\t".join(r) + "\n" for r in lists.get(mark, [])))
        paths.append(path)
    tf = tmp_path / "TF-List.txt"
    tf.write_text("".join("\t".join(r) + "\n" for r in tf_rows))
    return paths, tf


def test_na_file_counts_as_missing(tmp_path):
    """Treat a mark row whose file is NA as no experiment."""
    lists = {
        "DNase": [["ENCSR000AAA", "ENCFF000AAA", "cell_a"]],
        "H3K4me3": [
            ["ENCSR000AAB", "NA", "cell_a"],
            ["ENCSR000AAC", "NA", "cell_b"],
        ],
        "H3K27ac": [["ENCSR000AAD", "ENCFF000AAD", "cell_a"]],
        "CTCF": [["ENCSR000AAE", "ENCFF000AAE", "cell_a"]],
    }
    paths, tf = write_lists(tmp_path, lists)
    rows = mcl.master_rows(paths, tf)
    assert rows == [
        [
            "ENCSR000AAA",
            "ENCFF000AAA",
            "---",
            "---",
            "ENCSR000AAD",
            "ENCFF000AAD",
            "ENCSR000AAE",
            "ENCFF000AAE",
            "cell_a",
            "Group14",
            "no",
        ]
    ]
    assert mcl.output_name(rows[0]) == "ENCFF000AAA_ENCFF000AAD_ENCFF000AAE"
