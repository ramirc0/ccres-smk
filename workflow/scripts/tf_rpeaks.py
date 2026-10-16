"""Step 3b: reduce the TF peaks to representative peaks (rPeaks).

Port of 3_Curate-TF-Clusters.sh: peaks on chromosomes without '_' other than
chrEBV and chrM are sorted by position, then reduced by the merge, pick and
intersect loop on signalValue and name.
"""

import argparse

import _bed
import _registry
import polars as pl


def build_parser():
    """Build the tf_rpeaks.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tf_rpeaks.py flags.
    """
    p = argparse.ArgumentParser(description="Pick representative TF peaks.")
    p.add_argument("--peaks", required=True, help="Peak Parquet of tf_peaks.py.")
    p.add_argument("-o", "--output", required=True, help="rPeak Parquet to write.")
    return p


def main(argv=None):
    """Write the representative peaks in the order the loop appends them.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    peaks = (
        pl.read_parquet(args.peaks)
        .filter(_bed.primary(pl.col("chrom")))
        .sort(*_bed.COORDS, "name")
    )
    _bed.write_parquet(
        _registry.representative_peaks(peaks, "signal", "name"), args.output
    )


if __name__ == "__main__":
    main()
