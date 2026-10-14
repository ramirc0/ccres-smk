"""Step 2a: concatenate the processed DHSs and filter them.

Port of the first part of 2_Create-rDHSs.sh. `cat output.*` concatenated the
step 1 files in en_US glob order of their original names, output.<signal>.<row>.
The filter keeps DHSs at least 150 bp wide, with signal at least the cutoff
and FDR at most 0.001, on chromosomes without '_' other than chrEBV and chrM.
"""

import argparse
import locale
import shutil
from pathlib import Path

import _bed
import polars as pl

MIN_WIDTH = 150
MAX_FDR = 0.001
COLUMNS = ["chrom", "start", "end", "name", "fdr", "signal"]


def build_parser():
    """Build the combine_dhs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the combine_dhs.py flags.
    """
    p = argparse.ArgumentParser(description="Concatenate and filter DHSs.")
    p.add_argument("--hotspots", required=True, help="Hotspot list of the genome.")
    p.add_argument(
        "--processed", required=True, help="Directory of <signal>-<row>.bed files."
    )
    p.add_argument("--cutoff", type=float, required=True, help="Minimum signal.")
    p.add_argument("--all", required=True, help="Concatenated DHS BED to write.")
    p.add_argument("--filtered", required=True, help="Filtered DHS BED to write.")
    return p


def glob_order(hotspots):
    """Return the step 1 file stems in the original glob order.

    Parameters
    ----------
    hotspots : str
        Hotspot list. Column 3 is the signal file of each row.

    Returns
    -------
    list of str
        `<signal>-<row>` stems, ordered as en_US orders the original names,
        `output.<signal>.<row>`.
    """
    with open(hotspots) as f:
        pairs = [(line.split("\t")[2], row) for row, line in enumerate(f, start=1)]
    pairs.sort(key=lambda p: locale.strxfrm(f"output.{p[0]}.{p[1]}"))
    return [f"{signal}-{row}" for signal, row in pairs]


def main(argv=None):
    """Write the concatenated and the filtered DHSs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.all, "wb") as out:
        for stem in glob_order(args.hotspots):
            with open(Path(args.processed, f"{stem}.bed"), "rb") as f:
                shutil.copyfileobj(f, out, 1 << 24)

    dhs = _bed.read_tsv(args.all, COLUMNS)
    dhs = dhs.with_columns(
        signal_value=_bed.python_floats(dhs["signal"]),
        fdr_value=_bed.python_floats(dhs["fdr"]),
    )
    width = pl.col("end").cast(pl.Int64) - pl.col("start").cast(pl.Int64)
    keep = (
        _bed.primary(pl.col("chrom"))
        & (width >= MIN_WIDTH)
        & (pl.col("signal_value") >= args.cutoff)
        & (pl.col("fdr_value") <= MAX_FDR)
    )
    _bed.write_tsv(dhs.filter(keep).select(COLUMNS), args.filtered)


if __name__ == "__main__":
    main()
