"""Pieces shared by the cell type-agnostic (step 6) and specific (step 7) calls."""

import _bed
import polars as pl

HIGH = 1.64
CENTER = 200
FAR = 2000
TSS = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "transcript": pl.String,
    "score": pl.String,
    "strand": pl.String,
    "gene": pl.String,
}


def read_tss(path):
    """Read TSSs as the 1 bp intervals of their bases.

    TSS-Basic.bed writes a TSS's 1-based position p as `p p`, an empty
    interval. bedtools widened it to [p - 1, p + 1).

    Parameters
    ----------
    path : str
        TSS BED.

    Returns
    -------
    polars.DataFrame
        TSS rows with start p - 1 and end p.
    """
    return _bed.read_bed(path, TSS).with_columns(start=pl.col("end") - 1)


def high(mark):
    """Tell whether a z-score column is above 1.64.

    Parameters
    ----------
    mark : str
        Column of z-scores. Null for a mark without data.

    Returns
    -------
    polars.Expr
        True above 1.64, false otherwise, null included.
    """
    return (pl.col(mark) > HIGH).fill_null(False)


def label(rest):
    """Classify DNase-high regions by group and mark.

    A TSS region tests H3K4me3, then H3K27ac. Others test H3K27ac, then
    H3K4me3. CTCF comes last.

    Parameters
    ----------
    rest : polars.Expr
        Label of a region high in no mark.

    Returns
    -------
    polars.Expr
        PLS, pELS, dELS, CA-H3K4me3, CA-CTCF or `rest`, from the columns
        group, h3k4me3, h3k27ac and ctcf.
    """
    tss = pl.col("group") == "tss"
    return (
        pl.when(tss & high("h3k4me3"))
        .then(pl.lit("PLS"))
        .when(high("h3k27ac"))
        .then(
            pl.when(pl.col("group") == "distal")
            .then(pl.lit("dELS"))
            .otherwise(pl.lit("pELS"))
        )
        .when(~tss & high("h3k4me3"))
        .then(pl.lit("CA-H3K4me3"))
        .when(high("ctcf"))
        .then(pl.lit("CA-CTCF"))
        .otherwise(rest)
    )


def with_z(regions, key, tables):
    """Join z-score tables onto regions.

    Parameters
    ----------
    regions : polars.DataFrame
        Regions with a key column.
    key : str
        Column matched against each table's name.
    tables : dict of str to str or None
        Path of each mark's Parquet of name and z, or None for a mark without
        data.

    Returns
    -------
    polars.DataFrame
        The regions with one z column per mark, null without data.
    """
    for mark, path in tables.items():
        if path is None:
            regions = regions.with_columns(pl.lit(None, pl.Float64).alias(mark))
        else:
            z = pl.read_parquet(path, columns=["name", "z"]).rename(
                {"name": key, "z": mark}
            )
            regions = regions.join(z, on=key, how="left", maintain_order="left")
    return regions


def tss_groups(regions, tss, window):
    """Split regions into TSS, proximal and distal groups.

    Port of the shared block of 6_Classify-cCREs.sh and Split-cCREs.*.sh,
    with calculate-center-distance.py. TSS regions overlap a TSS. Proximal
    regions overlap a TSS +-2 kb window. Distal regions overlap none. A
    proximal region whose center lies within 200 bp of the nearest TSS start
    moves to TSS. One over 2 kb away moves to distal. The original measured
    from the first closest TSS bedtools printed, by edge distance.

    Parameters
    ----------
    regions, tss, window : polars.DataFrame
        Regions, TSSs and their windows, each with chrom, start and end.

    Returns
    -------
    polars.Series
        Group of each region: tss, prox or distal.
    """
    on_tss = _bed.overlaps_any(regions, tss)
    in_window = _bed.overlaps_any(regions, window)
    prox = regions.with_row_index("row").filter(~on_tss & in_window)
    starts = tss.select("chrom", tss=pl.col("start").cast(pl.Float64)).sort(
        "chrom", "tss"
    )
    nearest = (
        prox.select("row", "chrom", center=(pl.col("start") + pl.col("end")) / 2)
        .sort("chrom", "center")
        .join_asof(
            starts,
            left_on="center",
            right_on="tss",
            by="chrom",
            strategy="nearest",
            check_sortedness=False,
        )
        .select("row", distance=(pl.col("tss") - pl.col("center")).abs())
    )
    distance = regions.with_row_index("row").join(
        nearest, on="row", how="left", maintain_order="left"
    )["distance"]
    return pl.select(
        group=pl.when(on_tss | (distance <= CENTER))
        .then(pl.lit("tss"))
        .when(~in_window | (distance > FAR))
        .then(pl.lit("distal"))
        .otherwise(pl.lit("prox"))
    ).to_series()
