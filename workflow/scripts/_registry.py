"""Representative-peak selection (steps 2 and 3) and region accessioning."""

import _bed
import polars as pl

GENOME_ID = {"hg19": "EH37", "hg38": "EH38", "mm10": "EM10"}
MODE_ID = {"rDHS": "D", "cCRE": "E", "TF": "F"}


def representative_peaks(peaks, signal, name):
    """Pick representative peaks like the merge, pick and intersect loop.

    Each pass merges the remaining peaks (`bedtools merge`, book-ended
    included), picks the first peak of highest signal per cluster
    (`pick-best-peak.py`), keeps every remaining peak with a picked name
    (the awk join) and drops every remaining peak that overlaps a pick
    (`bedtools intersect -v`). Passes repeat until no peak remains. A peak
    overlaps only peaks of its own cluster, so the drop joins on the cluster.

    Parameters
    ----------
    peaks : polars.DataFrame
        Peaks with chrom, start and end, sorted by chrom and start.
    signal, name : str
        Columns of the signal, float, and of the peak name.

    Returns
    -------
    polars.DataFrame
        The representative peaks, pass by pass, in row order within a pass.
    """
    if peaks[signal].is_nan().any():
        raise ValueError("NaN signal: Python's max() would depend on its position.")
    remaining = peaks.with_row_index("_row")
    picks = []
    while not remaining.is_empty():
        remaining = remaining.with_columns(_cluster=_bed.clusters())
        best = remaining.group_by("_cluster").agg(
            pl.col(name).get(pl.col(signal).arg_max())
        )
        chosen = remaining.filter(pl.col(name).is_in(best[name].implode()))
        picks.append(chosen)
        hit = remaining.join(
            chosen.select("_cluster", pick_start="start", pick_end="end"),
            on="_cluster",
        ).filter(
            (pl.col("start") < pl.col("pick_end"))
            & (pl.col("end") > pl.col("pick_start"))
        )
        remaining = remaining.filter(~pl.col("_row").is_in(hit["_row"].implode()))
    return pl.concat(picks).drop("_row", "_cluster") if picks else peaks.clear()


def accessions(regions, genome, mode, previous=None):
    """Assign accessions to regions like make-region-accession.py.

    A region with the chrom, start and end of a previous region reuses that
    accession, the last previous row winning, with EH37 renamed EH38. The
    script keyed regions on `chrom + start + end`, so distinct coordinates
    could collide. Other regions get new numbers above the previous maximum,
    in row order.

    Parameters
    ----------
    regions : polars.DataFrame
        Regions with chrom, start and end, in output order.
    genome : str
        `hg19`, `hg38` or `mm10`.
    mode : str
        `rDHS`, `cCRE` or `TF`.
    previous : polars.DataFrame, optional
        Previous regions with chrom, start, end and accession.

    Returns
    -------
    polars.Series
        Accession of each region, named accession.
    """
    if previous is None:
        previous = pl.DataFrame(
            schema={**regions.select(_bed.COORDS).schema, "accession": pl.String}
        )
    top = (
        previous["accession"].str.split(MODE_ID[mode]).list.last().cast(pl.Int64).max()
    )
    known = previous.select(*_bed.COORDS, "accession").unique(
        _bed.COORDS, keep="last", maintain_order=True
    )
    found = regions.select(_bed.COORDS).join(
        known, on=_bed.COORDS, how="left", maintain_order="left"
    )["accession"]
    number = found.is_null().cum_sum() + (top or 0)
    new = GENOME_ID[genome] + MODE_ID[mode] + number.cast(pl.String).str.zfill(7)
    reused = found.str.replace("EH37", "EH38", literal=True)
    return pl.select(accession=pl.coalesce(reused, new)).to_series()
