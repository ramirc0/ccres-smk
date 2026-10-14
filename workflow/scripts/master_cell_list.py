"""Step 7a: the master cell list, one row per biosample with its mark experiments.

Port of match-biosamples.py and the TF flag of 7_Call-Cell-Type-Specific.sh.
Experiments join on their biosample key (Biosample-Keys.txt): the set of
biosamples their file came from, compared whole. The original joined them on
the biosample names in the lists, which spell one biosample several ways and
can name a pooled set after one donor. Rows come in order of first appearance
across the DNase, H3K4me3, H3K27ac and CTCF lists, and take the name of their
first experiment. Of two experiments of one mark under a key, the later
released wins, then the lower accession. A row whose file is NA counts as no
experiment. A TF experiment counts for the row whose key equals its own,
tagging modifications aside. Biosamples with TF data come first, flagged yes,
then the rest, flagged no.
"""

import argparse

MARKS = ("DNase", "H3K4me3", "H3K27ac", "CTCF")
GROUPS = {
    "YYYY": "Group1",
    "YYYN": "Group2",
    "YYNY": "Group3",
    "YYNN": "Group4",
    "YNYN": "Group5",
    "YNNY": "Group6",
    "YNNN": "Group7",
    "NYYN": "Group8",
    "NYNY": "Group9",
    "NYNN": "Group10",
    "NNYN": "Group11",
    "NNNY": "Group12",
    "NYYY": "Group13",
    "YNYY": "Group14",
    "NNYY": "Group15",
}
MISSING = "---"


def build_parser():
    """Build the master_cell_list.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the master_cell_list.py flags.
    """
    p = argparse.ArgumentParser(description="Build the master cell list.")
    for mark in MARKS:
        p.add_argument(f"--{mark.lower()}", required=True, help=f"{mark} list.")
    p.add_argument("--tf-list", required=True, help="Cell type TF experiment list.")
    p.add_argument("--keys", required=True, help="Biosample-Keys.txt.")
    p.add_argument("-o", "--output", required=True, help="Master cell list to write.")
    return p


def read_rows(path):
    """Read a tab-separated list as match-biosamples.py does.

    Parameters
    ----------
    path : str
        List file.

    Returns
    -------
    list of list of str
        Fields of each line after stripping trailing whitespace.
    """
    with open(path) as f:
        return [line.rstrip().split("\t") for line in f]


def read_keys(path):
    """Read Biosample-Keys.txt.

    Parameters
    ----------
    path : str
        Rows of experiment, file, key, key without tagging modifications, and
        release date.

    Returns
    -------
    dict of tuple of (str, str) to tuple of (str, str, str)
        The keys and release date of each (experiment, file).
    """
    return {(r[0], r[1]): (r[2], r[3], r[4]) for r in read_rows(path)}


def tf_experiments(tf_list, keys):
    """Group the TF list's experiments by biosample key.

    Parameters
    ----------
    tf_list : str
        TF experiment list. Columns 1-2 are experiment and peak file.
    keys : dict
        Output of `read_keys`.

    Returns
    -------
    dict of str to list of tuple of (str, str)
        (experiment, peak file) per key, tagging modifications aside, in list
        order.
    """
    grouped = {}
    for exp, peak, *_ in read_rows(tf_list):
        grouped.setdefault(keys[(exp, peak)][1], []).append((exp, peak))
    return grouped


def master_rows(mark_lists, tf_list, keys):
    """Return the master cell list rows, each with its biosample key.

    Parameters
    ----------
    mark_lists : list of str
        DNase, H3K4me3, H3K27ac and CTCF lists, in that order. Columns:
        experiment, file, biosample name.
    tf_list : str
        TF experiment list.
    keys : str
        Biosample-Keys.txt.

    Returns
    -------
    list of list of str
        Experiment and file of each mark (--- when absent), biosample name,
        group, TF flag, then the key, which the written list leaves out.
    """
    key_of = read_keys(keys)
    candidates, names = {}, {}
    for mode, path in enumerate(mark_lists):
        for exp, acc, cell, *_ in read_rows(path):
            if acc == "NA":
                continue
            key, _, released = key_of[(exp, acc)]
            names.setdefault(key, cell)
            marks = candidates.setdefault(key, [[] for _ in MARKS])
            marks[mode].append((exp, acc, released))
    with_tf = tf_experiments(tf_list, key_of)
    rows = []
    for key, marks in candidates.items():
        files = []
        for mark, found in zip(MARKS, marks):
            # max keeps the first of equal dates, so the lowest accession.
            best = max(sorted(found), key=lambda c: c[2]) if found else None
            if len(found) > 1:
                print(f"{mark} of {names[key]}: {best[0]} of {[c[0] for c in found]}")
            files += best[:2] if best else [MISSING, MISSING]
        group = GROUPS["".join("N" if f == MISSING else "Y" for f in files[::2])]
        rows.append([*files, names[key], group, "yes" if key in with_tf else "no", key])
    return [r for r in rows if r[10] == "yes"] + [r for r in rows if r[10] == "no"]


def output_name(row):
    """Return the name of a biosample's cCRE file, as the Split scripts build it.

    Parameters
    ----------
    row : list of str
        Master cell list row.

    Returns
    -------
    str
        File accessions joined by `_`, absent marks dropped. Without DNase,
        the name starts at H3K4me3.
    """
    if row[0] != MISSING:
        return "_".join(row[1:8:2]).replace("_" + MISSING, "")
    joined = "_".join(row[3:8:2])
    return joined.replace("_" + MISSING, "").replace(MISSING + "_", "")


def main(argv=None):
    """Write the master cell list.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    lists = [getattr(args, mark.lower()) for mark in MARKS]
    with open(args.output, "w") as f:
        f.writelines(
            "\t".join(r[:-1]) + "\n"
            for r in master_rows(lists, args.tf_list, args.keys)
        )


if __name__ == "__main__":
    main()
