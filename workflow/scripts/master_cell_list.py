"""Step 7a: the master cell list, one row per biosample with its mark experiments.

Port of match-biosamples.py and the TF flag of 7_Call-Cell-Type-Specific.sh.
Biosamples come in order of first appearance across the DNase, H3K4me3,
H3K27ac and CTCF lists. A later experiment of the same mark replaces an
earlier one. A row whose file is NA counts as no experiment. Biosamples with
TF data come first, flagged yes, then the rest, flagged no.
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


def master_rows(mark_lists, tf_list):
    """Return the master cell list rows.

    Parameters
    ----------
    mark_lists : list of str
        DNase, H3K4me3, H3K27ac and CTCF lists, in that order. Columns:
        experiment, file, biosample.
    tf_list : str
        TF experiment list. The last column is the biosample.

    Returns
    -------
    list of list of str
        Experiment and file of each mark (--- when absent), biosample, group
        and TF flag.
    """
    flags, files = {}, {}
    for mode, path in enumerate(mark_lists):
        for exp, acc, cell, *_ in read_rows(path):
            if acc == "NA":
                continue
            if cell not in flags:
                flags[cell] = ["N"] * len(MARKS)
                files[cell] = [MISSING] * (2 * len(MARKS))
            flags[cell][mode] = "Y"
            files[cell][2 * mode : 2 * mode + 2] = [exp, acc]
    with open(tf_list) as f:
        with_tf = {line.rstrip("\n").split("\t")[-1] for line in f}
    rows = [[*files[c], c, GROUPS["".join(flags[c])]] for c in flags]
    return [r + ["yes"] for r in rows if r[8] in with_tf] + [
        r + ["no"] for r in rows if r[8] not in with_tf
    ]


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
        f.writelines("\t".join(r) + "\n" for r in master_rows(lists, args.tf_list))


if __name__ == "__main__":
    main()
