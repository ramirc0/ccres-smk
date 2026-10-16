"""Step 3a: concatenate the TF peaks, resized to 150-350 bp, with their summits.

Port of the peak loop of 3_Curate-TF-Clusters.sh. Each experiment of the
filtered list adds its narrowPeak rows, in list order, with its experiment,
target and biosample, and the name `<experiment>-<line>`. The summit is start
+ column 10. A peak under 150 bp becomes summit +-75 and one over 350 bp
summit +-175. Other peaks pass unchanged.

The original also wrote tmp-summits, [start + column 10, +1) of each resized
row, with the old offset, which moved a resized peak's point off its summit.
tf_clusters takes the summit column instead.
"""

import argparse

import _bed
import _encode
import polars as pl

NARROWPEAK = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "name": pl.String,
    "score": pl.String,
    "strand": pl.String,
    "signal": pl.Float64,
    "p": pl.Float64,
    "q": pl.Float64,
    "offset": pl.Int64,
}
MIN_WIDTH, MAX_WIDTH = 150, 350


def build_parser():
    """Build the tf_peaks.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tf_peaks.py flags.
    """
    p = argparse.ArgumentParser(description="Concatenate resized TF peaks.")
    p.add_argument("--tf-list", required=True, help="Filtered TF experiment list.")
    p.add_argument("--encode-root", required=True, help="Local ENCODE mirror.")
    p.add_argument("--download-dir", required=True, help="Downloaded ENCODE files.")
    p.add_argument("-o", "--output", required=True, help="Peak Parquet to write.")
    return p


def resized(peaks, exp, target, biosample):
    """Resize one experiment's peaks and add the experiment columns.

    Parameters
    ----------
    peaks : polars.DataFrame
        The narrowPeak rows.
    exp, target, biosample : str
        Experiment, target and biosample of every row.

    Returns
    -------
    polars.DataFrame
        Chrom, start, end, summit, signal, experiment, target, biosample and
        name of each peak.
    """
    summit = pl.col("start") + pl.col("offset")
    width = pl.col("end") - pl.col("start")
    half = (
        pl.when(width < MIN_WIDTH)
        .then(MIN_WIDTH // 2)
        .when(width > MAX_WIDTH)
        .then(MAX_WIDTH // 2)
    )
    return peaks.select(
        "chrom",
        start=pl.coalesce(summit - half, "start"),
        end=pl.coalesce(summit + half, "end"),
        summit=summit,
        signal="signal",
        exp=pl.lit(exp),
        target=pl.lit(target),
        biosample=pl.lit(biosample),
        name=pl.format("{}-{}", pl.lit(exp), pl.int_range(1, pl.len() + 1)),
    )


def main(argv=None):
    """Write the peaks of every experiment of the list.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.tf_list) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    peaks = [
        resized(
            _bed.read_bed(
                _encode.encode_path(
                    args.encode_root, args.download_dir, exp, peak, "bed.gz"
                ),
                NARROWPEAK,
            ),
            exp,
            target,
            biosample,
        )
        for exp, peak, biosample, target, *_ in rows
    ]
    _bed.write_parquet(pl.concat(peaks), args.output)


if __name__ == "__main__":
    main()
