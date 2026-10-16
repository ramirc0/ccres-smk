"""Differential tests of _bed against bedtools 2.30.0."""

import _bed
import numpy as np
import polars as pl
import pytest

CHROMS = ["chr1", "chr2", "chr10", "chrX", "chrM", "chr1_KI270706v1_random"]
FDRS = [0.0468394, 1.82817e-05, 1e-300, 0.000999, 2.5e-07]


def random_bed(rng, n, span=3000, max_len=300):
    """Draw dense random intervals so that edge cases occur often.

    Parameters
    ----------
    rng : numpy.random.Generator
        Random source.
    n : int
        Number of intervals.
    span : int
        Starts are drawn from [0, span).
    max_len : int
        Lengths are drawn from [1, max_len).

    Returns
    -------
    polars.DataFrame
        Columns chrom, start, end, name, fdr.
    """
    start = rng.integers(0, span, n)
    return pl.DataFrame(
        {
            "chrom": rng.choice(CHROMS, n).tolist(),
            "start": start,
            "end": start + rng.integers(1, max_len, n),
            "name": [f"r{i}" for i in range(n)],
            "fdr": rng.choice(FDRS, n),
        }
    )


def write(df, tmp_path, name):
    """Write a table as BED text and return its path.

    Parameters
    ----------
    df : polars.DataFrame
        Table to write.
    tmp_path : pathlib.Path
        Directory.
    name : str
        File name.

    Returns
    -------
    pathlib.Path
        Written file.
    """
    path = tmp_path / name
    _bed.write_bed(df, path)
    return path


def text(df):
    """Render a table as BED text.

    Parameters
    ----------
    df : polars.DataFrame
        Table.

    Returns
    -------
    str
        Tab-separated lines, each ending in a newline.
    """
    return df.write_csv(separator="\t", include_header=False, quote_style="never")


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("distance", [0, 1])
def test_clusters(bedtools, tmp_path, seed, distance):
    """Match `bedtools merge -d` clusters, with collapse and min."""
    rng = np.random.default_rng(seed)
    df = random_bed(rng, 2000, span=20000, max_len=200).sort("chrom", "start")
    path = write(df, tmp_path, "a.bed")
    ours = (
        df.group_by(_bed.clusters(distance).alias("cluster"), maintain_order=True)
        .agg(
            pl.col("chrom").first(),
            pl.col("start").first(),
            pl.col("end").max(),
            pl.col("name").str.join(","),
            pl.col("fdr").min().map_elements(lambda v: f"{v:.10g}", pl.String),
        )
        .drop("cluster")
    )
    expected = bedtools(
        "merge", "-i", path, "-d", distance, "-c", "4,5", "-o", "collapse,min"
    )
    assert text(ours) == expected


@pytest.mark.parametrize("seed", range(4))
def test_overlap_pairs(bedtools, tmp_path, seed):
    """overlap_pairs gives the pairs of `bedtools intersect -wa -wb`."""
    rng = np.random.default_rng(seed)
    a = random_bed(rng, 1000)
    b = random_bed(rng, 1000, max_len=60)
    pa, pb = write(a, tmp_path, "a.bed"), write(b, tmp_path, "b.bed")
    pairs = _bed.overlap_pairs(a, b)
    b_hits = b[pairs["b"]].rename(lambda c: f"b_{c}")
    ours = text(pl.concat([a[pairs["a"]], b_hits], how="horizontal"))
    expected = bedtools("intersect", "-wa", "-wb", "-a", pa, "-b", pb)
    assert sorted(ours.splitlines()) == sorted(expected.splitlines())
    assert pairs.equals(pairs.sort("a", "b"))


@pytest.mark.parametrize("seed", range(4))
def test_overlaps_any(bedtools, tmp_path, seed):
    """overlaps_any reproduces `bedtools intersect -u` and `-v`."""
    rng = np.random.default_rng(seed)
    a = random_bed(rng, 1500)
    b = random_bed(rng, 1500, max_len=60)
    pa, pb = write(a, tmp_path, "a.bed"), write(b, tmp_path, "b.bed")
    hit = _bed.overlaps_any(a, b)
    assert text(a.filter(hit)) == bedtools("intersect", "-u", "-a", pa, "-b", pb)
    assert text(a.filter(~hit)) == bedtools("intersect", "-v", "-a", pa, "-b", pb)


def test_overlaps_handle_empty_tables():
    """An empty side gives no pairs."""
    a = random_bed(np.random.default_rng(0), 10)
    assert _bed.overlap_pairs(a, a.clear()).is_empty()
    assert not _bed.overlaps_any(a, a.clear()).any()
    assert _bed.overlaps_any(a.clear(), a).is_empty()


def test_overlaps_reject_zero_length():
    """Zero-length intervals raise. bedtools widened them to 2 bp."""
    a = pl.DataFrame({"chrom": ["chr1"], "start": [5], "end": [5]})
    with pytest.raises(ValueError):
        _bed.overlap_pairs(a, a)


def test_read_bed_types_columns(tmp_path):
    """read_bed casts each field and reads an empty file as zero rows."""
    path = tmp_path / "a.bed"
    path.write_text("chr1\t10\t20\t1.0e-05\nchr2\t5\t7\t0.0468394\n")
    schema = {"chrom": pl.String, "start": pl.Int64, "end": pl.Int64, "v": pl.Float64}
    df = _bed.read_bed(path, schema)
    assert df.schema == pl.Schema(schema)
    assert df.row(0) == ("chr1", 10, 20, 1e-05)
    empty = tmp_path / "empty.bed"
    empty.write_text("")
    assert _bed.read_bed(empty, schema).shape == (0, 4)


def test_primary_chromosomes_judge_the_chromosome_alone():
    """Drop chrM, chrEBV and names with '_', whatever the other fields hold."""
    rows = pl.DataFrame(
        {
            "chrom": ["chr1", "chrM", "chrEBV", "chr1_KI270706v1_random", "chr2"],
            "name": ["peak", "peak", "peak", "peak", "chrM_peak"],
        }
    )
    kept = rows.filter(_bed.primary(pl.col("chrom")))
    assert kept["chrom"].to_list() == ["chr1", "chr2"]
