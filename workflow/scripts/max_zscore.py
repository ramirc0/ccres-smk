"""Step 5: the highest z-score of each anchor across one mark's experiments.

Port of 5_Determine-Max-Zscores.sh with select-max-zscore.py. The files are
read in the en_US glob order of their names. The first seeds each anchor's
value, and a later value replaces it only if it is greater. Rows keep the
first file's order. Values print with Python's str().
"""

import argparse
import locale
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import _bed
import numpy as np


def build_parser():
    """Build the max_zscore.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the max_zscore.py flags.
    """
    p = argparse.ArgumentParser(description="Max z-score per anchor for one mark.")
    p.add_argument("--mark-list", required=True, help="Experiment list of the mark.")
    p.add_argument("--signal", required=True, help="Directory of <exp>-<file>.txt.")
    p.add_argument("-o", "--output", required=True, help="maxZ table to write.")
    p.add_argument("--threads", type=int, default=8, help="Files read at once.")
    return p


def read_z(path):
    """Read the names and z-scores of one signal file.

    Parameters
    ----------
    path : pathlib.Path
        Output of signal_zscore.py.

    Returns
    -------
    tuple of (polars.Series, numpy.ndarray)
        Names and float64 z-scores. Both empty for an empty file.
    """
    table = _bed.read_tsv(path, ["name", "z", "signal", "rank"])
    return table["name"], _bed.python_floats(table["z"])


def main(argv=None):
    """Write each anchor's maximum z-score.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.mark_list) as f:
        names = [
            "{}-{}.txt".format(*line.split("\t")[:2]) for line in f if line.strip()
        ]
    paths = [Path(args.signal, n) for n in sorted(names, key=locale.strxfrm)]

    with ThreadPoolExecutor(args.threads) as pool:
        tables = pool.map(read_z, paths)
        anchors, best = next(tables)
        if anchors.is_empty():
            raise ValueError(f"{paths[0]} is empty. The original would crash here.")
        for names_i, z in tables:
            if names_i.is_empty():
                continue
            if not names_i.equals(anchors):
                raise ValueError("Signal files list different anchors.")
            # A NaN never wins `>`, as in the original.
            best = np.where(z > best, z, best)

    with open(args.output, "w") as f:
        f.writelines(
            f"{name}\t{value!r}\n"
            for name, value in zip(anchors.to_list(), best.tolist())
        )


if __name__ == "__main__":
    main()
