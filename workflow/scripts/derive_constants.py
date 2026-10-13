"""Derive the step 0 and step 2 constants from the data.

The release hardcodes min_p and dhs_signal_cutoff, and the run keeps those
values. This reports what the data give for each, next to the release value:

- min_p: the largest FDR order in the hotspot files, read as call_dhs reads it.
- dhs_signal_cutoff: the 10th percentile of the DHS-All signal. The release
  took 10% of 10,000,000 randomly selected DHSs and kept no seed.
"""

import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

import _encode
import numpy as np
import polars as pl
from call_dhs import KEEP, fdr_order


def build_parser():
    """Build the derive_constants.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the derive_constants.py flags.
    """
    p = argparse.ArgumentParser(description="Derive min_p and the DHS cutoff.")
    p.add_argument("--hotspots", required=True, help="Hotspot list.")
    p.add_argument("--encode-root", required=True, help="Local ENCODE mirror.")
    p.add_argument("--download-dir", required=True, help="Downloaded ENCODE files.")
    p.add_argument("--dhs-all", required=True, help="DHS-All BED of step 2.")
    p.add_argument("--min-p", required=True, help="Release min_p.")
    p.add_argument("--cutoff", required=True, help="Release dhs_signal_cutoff.")
    p.add_argument("--threads", type=int, default=1, help="Files read at once.")
    p.add_argument("-o", "--output", required=True, help="TSV to write.")
    return p


def column(path, index):
    """Read one column of a tab-separated file as text.

    Parameters
    ----------
    path : str
        File to read. Gzip input is decompressed.
    index : int
        0-based column index.

    Returns
    -------
    polars.Series
        Column values as written.
    """
    df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        columns=[index],
        quote_char=None,
        infer_schema=False,
    )
    return df.to_series()


def max_order(path):
    """Return the largest FDR order of a hotspot file.

    Parameters
    ----------
    path : str
        Hotspot2 BED, gzipped. Column 5 holds the FDR.

    Returns
    -------
    int
        Largest order below `KEEP`, 0 if there is none.
    """
    orders = map(fdr_order, column(path, 4).unique().to_list())
    return max((order for order in orders if order != KEEP), default=0)


def main(argv=None):
    """Write each constant's release value next to its derived value.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.hotspots) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    paths = [
        _encode.encode_path(args.encode_root, args.download_dir, exp, peak, "bed.gz")
        for exp, peak, *_ in rows
    ]
    # Forking a process that runs polars' thread pool can deadlock.
    spawn = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(args.threads, mp_context=spawn) as pool:
        min_p = max(pool.map(max_order, paths))
    signal = column(args.dhs_all, 5).cast(pl.Float64).to_numpy()
    cutoff = np.percentile(signal, 10)
    with open(args.output, "w") as f:
        f.write("constant\trelease\tderived\tderivation\n")
        f.write(
            f"min_p\t{args.min_p}\t{min_p}\tlargest FDR order in the hotspot files\n"
        )
        f.write(
            f"dhs_signal_cutoff\t{args.cutoff}\t{cutoff:.6g}\t"
            "10th percentile of the DHS-All signal\n"
        )


if __name__ == "__main__":
    main()
