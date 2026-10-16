"""Step 2b: pick representative DHSs, accession them and keep the cDHS-supported ones.

Port of the rest of 2_Create-rDHSs.sh. The filtered DHSs are sorted by
position and reduced to representative DHSs (rDHSs) by the merge, pick and
intersect loop. The sorted rDHSs get accessions, reusing those of identical
earlier rDHSs. Each rDHS then falls in one cDHS class:

- complete: it covers a whole cDHS (`intersect -F 1 -u`).
- none: it overlaps no cDHS.
- 135: it overlaps a cDHS by at least 135 bp.
- less: the rest. It is kept if its signal is at least the 90th percentile of
  the none class's signal.

The summary holds the kept less, the 135 and the complete rDHSs, sorted by
position. The original printed the three classes one after another.
"""

import argparse

import _bed
import _registry
import polars as pl

ORDER = [*_bed.COORDS, "name"]
PREVIOUS = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "accession": pl.String,
    "version": pl.String,
}
CDHS = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    **{f"c{i}": pl.String for i in range(4, 12)},
}
MIN_OVERLAP = 135
PERCENTILE = 0.9


def build_parser():
    """Build the create_rdhs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the create_rdhs.py flags.
    """
    p = argparse.ArgumentParser(description="Create rDHSs from filtered DHSs.")
    p.add_argument("--filtered", required=True, help="Filtered DHS Parquet.")
    p.add_argument("--previous", required=True, help="Previous rDHSs, gzipped BED.")
    p.add_argument("--cdhs", required=True, help="Altius consensus DHS BED.")
    p.add_argument("--genome", required=True, help="hg38 or mm10.")
    p.add_argument("--prestam", required=True, help="Accessioned rDHSs to write.")
    p.add_argument("--complete", required=True, help="complete class to write.")
    p.add_argument("--none", required=True, help="none class to write.")
    p.add_argument("--summary", required=True, help="Kept rDHSs to write.")
    p.add_argument("--rdhs", required=True, help="rDHS BED to write.")
    return p


def cdhs_classes(rdhs, cdhs):
    """Measure each rDHS against the cDHSs.

    Parameters
    ----------
    rdhs, cdhs : polars.DataFrame
        Tables with chrom, start and end.

    Returns
    -------
    polars.DataFrame
        Per rDHS: whether it covers a whole cDHS, whether it overlaps one, and
        its largest overlap in bp.
    """
    pairs = _bed.overlap_pairs(rdhs, cdhs)
    a = rdhs.select("start", "end")[pairs["a"]]
    b = cdhs.select(b_start="start", b_end="end")[pairs["b"]]
    measured = (
        pl.concat([pairs.select(row="a"), a, b], how="horizontal")
        .group_by("row")
        .agg(
            covers=(
                (pl.col("b_start") >= pl.col("start"))
                & (pl.col("b_end") <= pl.col("end"))
            ).any(),
            overlap=(
                pl.min_horizontal("end", "b_end")
                - pl.max_horizontal("start", "b_start")
            ).max(),
        )
    )
    return (
        rdhs.with_row_index("row")
        .join(measured, on="row", how="left", maintain_order="left")
        .select(
            covers=pl.col("covers").fill_null(False),
            overlaps=pl.col("overlap").is_not_null(),
            overlap=pl.col("overlap").fill_null(0),
        )
    )


def main(argv=None):
    """Write the accessioned rDHSs, their cDHS classes and the kept rDHSs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    dhs = pl.read_parquet(args.filtered).sort(ORDER)
    rdhs = _registry.representative_peaks(dhs, "signal", "name").sort(ORDER)
    previous = _bed.read_bed(args.previous, PREVIOUS)
    zz = rdhs.with_columns(_registry.accessions(rdhs, args.genome, "rDHS", previous))
    _bed.write_parquet(zz, args.prestam)

    zz = pl.concat(
        [zz, cdhs_classes(zz, _bed.read_bed(args.cdhs, CDHS))], how="horizontal"
    )
    complete = zz.filter("covers")
    _bed.write_parquet(complete.select(*rdhs.columns, "accession"), args.complete)
    rest = zz.filter(~pl.col("covers"))
    none = rest.filter(~pl.col("overlaps"))
    _bed.write_parquet(none.select(*rdhs.columns, "accession"), args.none)

    wide = rest.filter(pl.col("overlap") >= MIN_OVERLAP)
    cutoff = none["signal"].quantile(PERCENTILE, interpolation="linear")
    less = rest.filter(
        pl.col("overlaps")
        & ~pl.col("name").is_in(wide["name"].implode())
        & (pl.col("signal") >= cutoff)
    )
    summary = pl.concat([less, wide, complete]).select(*rdhs.columns, "accession")
    summary = summary.sort(ORDER)
    _bed.write_parquet(summary, args.summary)
    _bed.write_bed(summary.select(*_bed.COORDS, "accession"), args.rdhs)


if __name__ == "__main__":
    main()
