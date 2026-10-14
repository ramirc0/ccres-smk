"""Step 3a: concatenate the TF peaks, resized to 150-350 bp, and their summits.

Port of the peak loop of 3_Curate-TF-Clusters.sh. Each experiment of the
filtered list adds its narrowPeak rows, in list order, with four more columns:
experiment, target, biosample and `<experiment>-<line>`. A peak under 150 bp
becomes summit +-75 and one over 350 bp summit +-175, where the summit is
start + column 10. Other peaks pass unchanged.

tmp-summits holds [start + column 10, +1) of each row. A resized peak's
column 10 becomes 75 or 175, so that point stays on its summit. The original
kept the old offset, which moved the point off the summit.
"""

import argparse

import _bed
import _encode
import polars as pl

NARROWPEAK = [f"c{i}" for i in range(1, 11)]
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
    p.add_argument("--bed", required=True, help="tmp.bed to write.")
    p.add_argument("--summits", required=True, help="tmp-summits to write.")
    return p


def resized(peaks, exp, target, biosample):
    """Resize one experiment's peaks and add the experiment columns.

    Parameters
    ----------
    peaks : polars.DataFrame
        The narrowPeak string columns c1-c10.
    exp, target, biosample : str
        Columns 11-13 of every row.

    Returns
    -------
    polars.DataFrame
        14 string columns, as the original awk printed them, except column
        10 of resized peaks.
    """
    if any((peaks[c] == "").any() for c in NARROWPEAK):
        raise ValueError(f"{exp}: empty field. awk would shift the columns.")
    start, end = pl.col("c2").cast(pl.Int64), pl.col("c3").cast(pl.Int64)
    summit = start + pl.col("c10").cast(pl.Int64)
    width = end - start
    half = (
        pl.when(width < MIN_WIDTH)
        .then(MIN_WIDTH // 2)
        .when(width > MAX_WIDTH)
        .then(MAX_WIDTH // 2)
    )
    # Unresized rows keep their text: awk printed them as $0.
    return peaks.with_columns(
        c2=pl.coalesce((summit - half).cast(pl.String), pl.col("c2")),
        c3=pl.coalesce((summit + half).cast(pl.String), pl.col("c3")),
        c10=pl.coalesce(half.cast(pl.String), pl.col("c10")),
        c11=pl.lit(exp),
        c12=pl.lit(target),
        c13=pl.lit(biosample),
        c14=pl.format("{}-{}", pl.lit(exp), pl.int_range(1, pl.len() + 1)),
    )


def summit_points(bed):
    """Return the tmp-summits rows of tmp.bed rows.

    Parameters
    ----------
    bed : polars.DataFrame
        14 string columns of tmp.bed.

    Returns
    -------
    polars.DataFrame
        Rows with columns 2-3 replaced by [c2 + c10, c2 + c10 + 1).
    """
    point = pl.col("c2").cast(pl.Int64) + pl.col("c10").cast(pl.Int64)
    return bed.with_columns(c2=point.cast(pl.String), c3=(point + 1).cast(pl.String))


def main(argv=None):
    """Write tmp.bed and tmp-summits for every experiment of the list.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.tf_list) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    with open(args.bed, "w") as bed_out, open(args.summits, "w") as summit_out:
        for exp, peak, biosample, target, *_ in rows:
            path = _encode.encode_path(
                args.encode_root, args.download_dir, exp, peak, "bed.gz"
            )
            bed = resized(_bed.read_tsv(path, NARROWPEAK), exp, target, biosample)
            for table, out in ((bed, bed_out), (summit_points(bed), summit_out)):
                out.writelines(f"{line}\n" for line in _bed.lines(table).to_list())


if __name__ == "__main__":
    main()
