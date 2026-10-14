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


def write_keys(tmp_path, rows):
    """Write a Biosample-Keys.txt.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to write in.
    rows : list of tuple of str
        Experiment, file, key, key without tagging modifications, release date.

    Returns
    -------
    pathlib.Path
        The file.
    """
    path = tmp_path / "Biosample-Keys.txt"
    path.write_text("".join("\t".join(r) + "\n" for r in rows))
    return path


def test_rows_join_on_biosample_key(tmp_path):
    """Join experiments whose keys match, whatever the names in the lists."""
    lists = {
        "DNase": [["ENCSR000AAA", "ENCFF000AAA", "cell_a"]],
        "H3K4me3": [["ENCSR000AAB", "ENCFF000AAB", "cell_a_spelled_otherwise"]],
        "H3K27ac": [["ENCSR000AAC", "ENCFF000AAC", "cell_a"]],
    }
    tf_rows = [["ENCSR000AAD", "ENCFF000AAD", "x", "CTCF", "0.01", "lab", "line", "NA"]]
    paths, tf = write_lists(tmp_path, lists, tf_rows)
    keys = write_keys(
        tmp_path,
        [
            ("ENCSR000AAA", "ENCFF000AAA", "k1", "k1", "2020-01-01"),
            ("ENCSR000AAB", "ENCFF000AAB", "k1", "k1", "2020-01-01"),
            ("ENCSR000AAC", "ENCFF000AAC", "k2", "k2", "2020-01-01"),
            ("ENCSR000AAD", "ENCFF000AAD", "k2_tagged", "k2", "2020-01-01"),
        ],
    )
    missing = ["---", "---"]
    assert mcl.master_rows(paths, tf, keys) == [
        [*missing * 2, "ENCSR000AAC", "ENCFF000AAC", *missing]
        + ["cell_a", "Group11", "yes", "k2"],
        ["ENCSR000AAA", "ENCFF000AAA", "ENCSR000AAB", "ENCFF000AAB", *missing * 2]
        + ["cell_a", "Group4", "no", "k1"],
    ]


def test_tf_flag_follows_biosample_key(tmp_path):
    """Flag a row when a TF experiment that passes step 3 shares its key."""
    lists = {
        "DNase": [
            ["ENCSR000AAA", "ENCFF000AAA", "cell_a"],
            ["ENCSR000AAB", "ENCFF000AAB", "cell_b"],
            ["ENCSR000AAC", "ENCFF000AAC", "cell_c"],
        ]
    }
    tf_rows = [
        [f"ENCSR000TF{i}", f"ENCFF000TF{i}", "x", "CTCF", frip, "lab", "line", "NA"]
        for i, frip in enumerate(["0.01", "0.001", "NA"])
    ]
    paths, tf = write_lists(tmp_path, lists, tf_rows)
    keys = write_keys(
        tmp_path,
        [
            *(
                (f"ENCSR000AA{c}", f"ENCFF000AA{c}", f"k{c}", f"k{c}", "2020-01-01")
                for c in "ABC"
            ),
            *(
                (f"ENCSR000TF{i}", f"ENCFF000TF{i}", "t", f"k{c}", "2020-01-01")
                for i, c in enumerate("ABC")
            ),
        ],
    )
    flags = {r[8]: r[10] for r in mcl.master_rows(paths, tf, keys)}
    assert flags == {"cell_a": "yes", "cell_b": "no", "cell_c": "no"}


def test_one_mark_twice_under_a_key(tmp_path):
    """Keep the later released experiment of a mark, then the lower accession."""
    lists = {
        "DNase": [
            ["ENCSR000AAA", "ENCFF000AAA", "cell_a"],
            ["ENCSR000AAB", "ENCFF000AAB", "cell_a"],
        ],
        "H3K4me3": [
            ["ENCSR000AAD", "ENCFF000AAD", "cell_a"],
            ["ENCSR000AAC", "ENCFF000AAC", "cell_a"],
        ],
    }
    paths, tf = write_lists(tmp_path, lists)
    keys = write_keys(
        tmp_path,
        [
            ("ENCSR000AAA", "ENCFF000AAA", "k", "k", "2016-01-01"),
            ("ENCSR000AAB", "ENCFF000AAB", "k", "k", "2012-01-01"),
            ("ENCSR000AAD", "ENCFF000AAD", "k", "k", "2020-10-02"),
            ("ENCSR000AAC", "ENCFF000AAC", "k", "k", "2020-10-02"),
        ],
    )
    (row,) = mcl.master_rows(paths, tf, keys)
    assert row[:4] == ["ENCSR000AAA", "ENCFF000AAA", "ENCSR000AAC", "ENCFF000AAC"]


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
    keys = write_keys(
        tmp_path,
        [
            (exp, f"ENCFF{exp[5:]}", "k", "k", "2020-01-01")
            for exp in ("ENCSR000AAA", "ENCSR000AAD", "ENCSR000AAE")
        ],
    )
    rows = mcl.master_rows(paths, tf, keys)
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
            "k",
        ]
    ]


def test_output_name_is_the_first_experiment():
    """Name a row's files after its first experiment: DNase, then the marks."""
    missing = ["---", "---"]
    with_dnase = ["ENCSR000AAA", "ENCFF000AAA", *missing, "ENCSR000AAB", "ENCFF000AAB"]
    without = [*missing * 2, "ENCSR000AAC", "ENCFF000AAC", "ENCSR000AAD", "ENCFF000AAD"]
    tail = ["cell", "Group5", "no", "key"]
    assert mcl.output_name([*with_dnase, *missing, *tail]) == "ENCSR000AAA"
    assert mcl.output_name([*without, *tail]) == "ENCSR000AAC"
