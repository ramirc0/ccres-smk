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

COLUMNS = [f"c{i}" for i in range(1, 15)]
MIN_EXPERIMENTS = 5


def build_parser():
    """Build the tf_clusters.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tf_clusters.py flags.
    """
    p = argparse.ArgumentParser(description="Curate TF clusters.")
    p.add_argument("--rpeaks", required=True, help="tmp.rPeaks.")
    p.add_argument("--summits", required=True, help="tmp.summits.")
    p.add_argument("--rdhs", required=True, help="rDHS BED of the genome.")
    p.add_argument("--genome", required=True, help="hg38 or mm10.")
    p.add_argument("--clusters", required=True, help="TF-Clusters.txt to write.")
    p.add_argument("--no-overlap", required=True, help="tmp.no-overlap to write.")
    p.add_argument("--summary", required=True, help="rTFBS-Summary.txt to write.")
    p.add_argument("--rtfbs", required=True, help="rTFBS.bed to write.")
    return p


def intervals(df, codes, start="c2", end="c3"):
    """Return the interval arrays of a table.

    Parameters
    ----------
    df : polars.DataFrame
        String columns.
    codes : numpy.ndarray
        Chromosome code of each row.
    start, end : str
        Coordinate columns.

    Returns
    -------
    tuple of numpy.ndarray
        Codes, starts and ends.
    """
    return codes, df[start].cast(pl.Int64).to_numpy(), df[end].cast(pl.Int64).to_numpy()


def clusters(rpeaks, summits):
    """Collect the summits of each rPeak and keep those with >= 5 experiments.

    Parameters
    ----------
    rpeaks, summits : polars.DataFrame
        14 string columns each.

    Returns
    -------
    polars.DataFrame
        Chrom, start, end, name, targets, biosamples, experiments and the
        experiment count, in rPeak order.
    """
    rp_codes, su_codes = _bed.shared_codes(rpeaks["c1"], summits["c1"])
    ia, ib, _ = _bed.overlap_pairs(
        intervals(rpeaks, rp_codes), intervals(summits, su_codes)
    )
    hits = pl.DataFrame(
        {
            "row": ia,
            "exp": summits["c11"].gather(ib),
            "target": summits["c12"].gather(ib),
            "cell": summits["c13"].gather(ib),
        }
    )
    joined = (
        hits.group_by("row", maintain_order=True)
        .agg(
            *(
                pl.col(c).unique().sort().str.join(";").alias(c)
                for c in ("target", "cell", "exp")
            ),
            count=pl.col("exp").n_unique(),
        )
        .filter(pl.col("count") >= MIN_EXPERIMENTS)
    )
    rows = joined["row"].to_numpy()
    return pl.DataFrame(
        {
            "chrom": rpeaks["c1"].gather(rows),
            "start": rpeaks["c2"].gather(rows),
            "end": rpeaks["c3"].gather(rows),
            "name": rpeaks["c14"].gather(rows),
            "target": joined["target"],
            "cell": joined["cell"],
            "exp": joined["exp"],
            "count": joined["count"].cast(pl.String),
        }
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
        _bed.read_tsv(args.rpeaks, COLUMNS), _bed.read_tsv(args.summits, COLUMNS)
    )
    _bed.write_tsv(tf, args.clusters)

    rdhs = _bed.read_tsv(args.rdhs, ["chrom", "start", "end", "accession"])
    tf_codes, rdhs_codes = _bed.shared_codes(tf["chrom"], rdhs["chrom"])
    off = ~_bed.overlaps_any(
        intervals(tf, tf_codes, "start", "end"),
        intervals(rdhs, rdhs_codes, "start", "end"),
    )
    lone = tf.filter(off)
    _bed.write_tsv(lone, args.no_overlap)

    regions = lone.select("chrom", "start", "end").iter_rows()
    summary = lone.with_columns(
        accession=pl.Series(
            _registry.accessions(regions, args.genome, "TF"), dtype=pl.String
        )
    )
    _bed.write_tsv(summary, args.summary)
    _bed.write_tsv(summary.select("chrom", "start", "end", "accession"), args.rtfbs)


if __name__ == "__main__":
    main()
