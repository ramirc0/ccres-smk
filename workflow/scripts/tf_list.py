r"""Step 3: keep the TF ChIP-seq experiments with FRiP > 0.003 and a peak file.

Ports `awk -F "\t" '{if ($5 > 0.003 && $2 != "NA") print $0}'`. awk compares
a FRiP that does not look numeric, such as NA, as a string, so it passes.
"""

import argparse
import re

# A field awk treats as a number.
AWK_NUMBER = re.compile(r"[ \t\n]*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?[ \t\n]*")


def build_parser():
    """Build the tf_list.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tf_list.py flags.
    """
    p = argparse.ArgumentParser(description="Filter the TF experiment list.")
    p.add_argument("-i", "--input", required=True, help="TF-List.All.txt.")
    p.add_argument("-o", "--output", required=True, help="Filtered list to write.")
    p.add_argument(
        "--min-frip", default="0.003", help="FRiP threshold, as written in awk."
    )
    return p


def awk_greater(field, constant):
    """Compare a field with a numeric constant as awk does.

    Parameters
    ----------
    field : str
        Field text.
    constant : str
        Numeric constant, as written in the awk program. awk prints it with
        CONVFMT (%.6g) when the comparison is a string comparison.

    Returns
    -------
    bool
        `field > constant` under awk semantics.
    """
    if AWK_NUMBER.fullmatch(field):
        return float(field) > float(constant)
    return field > f"{float(constant):.6g}"


def main(argv=None):
    """Write the experiments that pass the FRiP and peak-file filter.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.input) as f, open(args.output, "w") as o:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if awk_greater(fields[4], args.min_frip) and fields[1] != "NA":
                o.write(line)


if __name__ == "__main__":
    main()
