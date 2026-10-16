"""Typed BED tables: text and Parquet IO, and interval ops.

A table is a polars DataFrame with chrom (String), start and end (Int64)
columns, 0-based and half-open. Strings sort by bytes, as `LC_ALL=C sort`
does. Overlaps run on polars-bio. Zero-length intervals raise: bedtools widened
them to 2 bp, and the inputs that hold them read as 1 bp intervals instead.
"""

import subprocess
import tempfile
from pathlib import Path

import polars as pl
import polars_bio as pb

COORDS = ["chrom", "start", "end"]
AVERAGES = {
    "name": pl.String,
    "size": pl.Int64,
    "covered": pl.Int64,
    "sum": pl.Float64,
    "mean0": pl.Float64,
    "mean": pl.Float64,
}
# Fixed row groups make a table's Parquet bytes independent of how polars
# splits the work.
ROW_GROUP = 1 << 20


def read_bed(path, schema):
    """Read a headerless tab-separated file into typed columns.

    Parameters
    ----------
    path : str
        File to read. Gzip input is decompressed.
    schema : dict of str to polars.DataType
        Name and type of each column, in file order. Their number must match
        the file.

    Returns
    -------
    polars.DataFrame
        One column per field. An empty file gives zero rows.
    """
    df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        quote_char=None,
        infer_schema=False,
        empty_string_is_null=False,
        raise_if_empty=False,
    )
    if df.width == 0:
        return pl.DataFrame(schema=schema)
    if df.width != len(schema):
        raise ValueError(f"{path}: {df.width} columns, expected {len(schema)}.")
    df = df.rename(dict(zip(df.columns, schema)))
    # read_csv returns thousands of chunks per column, which make gathers slow.
    return df.cast(schema, strict=True).rechunk()


def write_bed(df, path):
    """Write a table as headerless tab-separated text.

    Parameters
    ----------
    df : polars.DataFrame
        Columns written in order. Floats print as the shortest text that
        reads back to the same double.
    path : str
        File to write.
    """
    df.write_csv(path, separator="\t", include_header=False, quote_style="never")


def write_parquet(df, path):
    """Write a table as zstd Parquet.

    Parameters
    ----------
    df : polars.DataFrame
        Table to write.
    path : str
        File to write.
    """
    df.write_parquet(path, compression="zstd", row_group_size=ROW_GROUP)


def bigwig_averages(regions, bigwig):
    """Average a bigWig over named regions with bigWigAverageOverBed.

    Parameters
    ----------
    regions : polars.DataFrame
        Chrom, start, end and name.
    bigwig : str
        Signal bigWig.

    Returns
    -------
    polars.DataFrame
        The tool's columns, in region order: name, size, covered, sum, mean0
        and mean. Its mean0 prints 6 significant digits.
    """
    with tempfile.TemporaryDirectory() as tmp:
        bed, tab = Path(tmp, "regions.bed"), Path(tmp, "out.tab")
        write_bed(regions.select(*COORDS, "name"), bed)
        subprocess.run(["bigWigAverageOverBed", bigwig, bed, tab], check=True)
        return read_bed(tab, AVERAGES)


def primary(chrom):
    """Tell which chromosomes steps 2 and 3 keep.

    The scripts dropped every line that contained chrEBV or chrM. Here only
    the chromosome counts.

    Parameters
    ----------
    chrom : polars.Expr
        Chromosome names.

    Returns
    -------
    polars.Expr
        False for chrEBV, chrM and names with '_'.
    """
    return ~chrom.str.contains("_", literal=True) & ~chrom.is_in(["chrEBV", "chrM"])


def clusters(distance=0):
    """Assign rows to the clusters of `bedtools merge -d distance`.

    A row joins the open cluster if its start is at most the furthest end of
    the cluster's rows plus `distance`, so book-ended rows merge at 0.

    Parameters
    ----------
    distance : int
        Largest gap that still merges.

    Returns
    -------
    polars.Expr
        Cluster of each row, from 1 in row order. Rows MUST be sorted by
        chrom and start, and each chrom's rows MUST be contiguous.
    """
    reach = pl.col("end").cum_max().shift().over("chrom")
    return (pl.col("start") > reach + distance).fill_null(True).cum_sum()


def _intervals(df):
    """Mark a table's intervals for polars-bio, with a row index.

    Parameters
    ----------
    df : polars.DataFrame
        Table with chrom, start and end.

    Returns
    -------
    polars.DataFrame
        Chrom, start, end and row, flagged 0-based.
    """
    if (df["end"] <= df["start"]).any():
        raise ValueError("Empty or zero-length interval.")
    out = df.select(*COORDS).with_row_index("row")
    out.config_meta.set(coordinate_system_zero_based=True)
    return out


def overlap_pairs(a, b):
    """List the overlapping rows of two tables, as `bedtools intersect -wa -wb`.

    Intervals overlap when they share a base. Book-ended ones do not.

    Parameters
    ----------
    a, b : polars.DataFrame
        Tables with chrom, start and end.

    Returns
    -------
    polars.DataFrame
        Row index in `a` and in `b` of each pair, sorted by both.
    """
    if a.is_empty() or b.is_empty():
        return pl.DataFrame(schema={"a": pl.UInt32, "b": pl.UInt32})
    pairs = pb.overlap(_intervals(a), _intervals(b), output_type="polars.DataFrame")
    return pairs.select(a="row_1", b="row_2").sort("a", "b")


def overlaps_any(a, b):
    """Flag rows of `a` that overlap a row of `b`, as `bedtools intersect -u`.

    Parameters
    ----------
    a, b : polars.DataFrame
        Tables with chrom, start and end.

    Returns
    -------
    polars.Series
        Flag per row of `a`.
    """
    hits = overlap_pairs(a, b)["a"]
    return pl.int_range(a.height, dtype=pl.UInt32, eager=True).is_in(hits.implode())
