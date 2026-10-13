"""Step 3b: reduce the TF peaks to representative peaks (rPeaks).

Port of 3_Curate-TF-Clusters.sh: peaks on chromosomes without '_' and on lines
without chrEBV or chrM are sorted with `sort -k1,1 -k2,2n`, then reduced by
the merge, pick and intersect loop on signalValue (column 7) and name
(column 14).
"""

import argparse

import _bed
import _registry
import polars as pl

COLUMNS = [f"c{i}" for i in range(1, 15)]


def build_parser():
    """Build the tf_rpeaks.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tf_rpeaks.py flags.
    """
    p = argparse.ArgumentParser(description="Pick representative TF peaks.")
    p.add_argument("--bed", required=True, help="tmp.bed of tf_peaks.py.")
    p.add_argument("-o", "--output", required=True, help="tmp.rPeaks to write.")
    return p


def main(argv=None):
    """Write the representative peaks in the order the loop appends them.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    peaks = _bed.read_tsv(args.bed, COLUMNS)
    line = pl.concat_str(COLUMNS, separator="\t")
    peaks = peaks.filter(
        ~pl.col("c1").str.contains("_", literal=True)
        & ~line.str.contains("chrEBV", literal=True)
        & ~line.str.contains("chrM", literal=True)
    )
    peaks = peaks[_bed.gnu_sort(peaks, text_keys=["c1"], numeric_keys=["c2"])]
    picks = _registry.representative_peaks(
        _bed.chrom_codes(peaks["c1"]),
        peaks["c2"].cast(pl.Int64).to_numpy(),
        peaks["c3"].cast(pl.Int64).to_numpy(),
        _bed.python_floats(peaks["c7"]),
        peaks["c14"].cast(pl.Categorical).to_physical().to_numpy(),
    )
    _bed.write_tsv(peaks[picks], args.output)


if __name__ == "__main__":
    main()
