"""Step 3c: keep rPeaks bound in at least 5 experiments, then accession those off rDHSs.

Port of the end of 3_Curate-TF-Clusters.sh with filter-tf-rpeaks.py. Each
rPeak collects the experiments, targets and biosamples of the summits it
overlaps (`bedtools intersect -wo`). rPeaks with summits from at least 5
experiments are the TF clusters, in rPeak order. Five experiments of one TF
count as five, though the paper's SI speaks of five transcription factors.
Clusters that overlap no rDHS get TF accessions in order.

The original joined each set in Python set() order, which varies between runs.
Here each set is joined in sorted order.
"""

import argparse

import _bed
import _registry
import polars as pl

MIN_EXPERIMENTS = 5
RDHS = {"chrom": pl.String, "start": pl.Int64, "end": pl.Int64, "accession": pl.String}


def build_parser():
    """Build the tf_clusters.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tf_clusters.py flags.
    """
    p = argparse.ArgumentParser(description="Curate TF clusters.")
    p.add_argument("--rpeaks", required=True, help="rPeak Parquet.")
    p.add_argument("--peaks", required=True, help="Peak Parquet, with summits.")
    p.add_argument("--rdhs", required=True, help="rDHS BED of the genome.")
    p.add_argument("--genome", required=True, help="hg38 or mm10.")
    p.add_argument("--clusters", required=True, help="TF clusters to write.")
    p.add_argument("--no-overlap", required=True, help="Clusters off rDHSs to write.")
    p.add_argument("--summary", required=True, help="Accessioned clusters to write.")
    p.add_argument("--rtfbs", required=True, help="rTFBS BED to write.")
    return p


def clusters(rpeaks, peaks):
    """Collect the summits of each rPeak and keep those with >= 5 experiments.

    Parameters
    ----------
    rpeaks : polars.DataFrame
        Representative peaks.
    peaks : polars.DataFrame
        All peaks, with summit, exp, target and biosample.

    Returns
    -------
    polars.DataFrame
        Chrom, start, end, name, targets, biosamples, experiments and the
        experiment count, in rPeak order.
    """
    summits = peaks.select(
        "chrom",
        "exp",
        "target",
        "biosample",
        start="summit",
        end=pl.col("summit") + 1,
    )
    pairs = _bed.overlap_pairs(rpeaks, summits)
    hits = summits.select("exp", "target", "biosample")[pairs["b"]].with_columns(
        row=pairs["a"]
    )
    joined = (
        hits.group_by("row")
        .agg(
            *(
                pl.col(c).unique().sort().str.join(";")
                for c in ("target", "biosample", "exp")
            ),
            count=pl.col("exp").n_unique(),
        )
        .filter(pl.col("count") >= MIN_EXPERIMENTS)
        .sort("row")
    )
    return pl.concat(
        [
            rpeaks.select(*_bed.COORDS, "name")[joined["row"]],
            joined.drop("row"),
        ],
        how="horizontal",
    )


def main(argv=None):
    """Write the TF clusters, those off rDHSs, and their accessioned summary.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    tf = clusters(
        pl.read_parquet(args.rpeaks),
        pl.read_parquet(
            args.peaks, columns=["chrom", "summit", "exp", "target", "biosample"]
        ),
    )
    _bed.write_parquet(tf, args.clusters)

    lone = tf.filter(~_bed.overlaps_any(tf, _bed.read_bed(args.rdhs, RDHS)))
    _bed.write_parquet(lone, args.no_overlap)

    summary = lone.with_columns(_registry.accessions(lone, args.genome, "TF"))
    _bed.write_parquet(summary, args.summary)
    _bed.write_bed(summary.select(*_bed.COORDS, "accession"), args.rtfbs)


if __name__ == "__main__":
    main()
