r"""Step 3: keep the TF ChIP-seq experiments with FRiP > 0.003 and a peak file.

The original filtered with `awk -F "\t" '{if ($5 > 0.003 && $2 != "NA") ...}'`.
awk compared a FRiP that does not look numeric, such as NA, as text, so it
passed. Here a FRiP must be a finite number.
"""

import argparse
import math

MIN_FRIP = 0.003


def build_parser():
    """Build the tf_list.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tf_list.py flags.
    """
    p = argparse.ArgumentParser(description="Filter the TF experiment list.")
    p.add_argument("-i", "--input", required=True, help="TF-List-All.txt.")
    p.add_argument("-o", "--output", required=True, help="Filtered list to write.")
    p.add_argument("--min-frip", type=float, default=MIN_FRIP, help="FRiP threshold.")
    return p


def passes(fields, min_frip=MIN_FRIP):
    """Return whether an experiment passes the FRiP and peak-file filter.

    Parameters
    ----------
    fields : list of str
        Fields of one TF-List-All.txt line.
    min_frip : float
        FRiP threshold.

    Returns
    -------
    bool
        True when FRiP is a finite number above the threshold and the
        experiment has a peak file.
    """
    try:
        frip = float(fields[4])
    except ValueError:
        return False
    return math.isfinite(frip) and frip > min_frip and fields[1] != "NA"


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
            if passes(line.rstrip("\n").split("\t"), args.min_frip):
                o.write(line)


if __name__ == "__main__":
    main()
