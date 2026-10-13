"""Differential tests of _bed against bedtools 2.30.0 and GNU sort under en_US."""

import _bed
import numpy as np
import polars as pl
import pytest

CHROMS = ["chr1", "chr2", "chr10", "chrX", "chrM", "chr1_KI270706v1_random"]
FDRS = ["0.0468394", "1.82817e-05", "3.2e-4944", "1e-300", "0.000999", "2.5e-07"]


def random_bed(rng, n, span=3000, max_len=300, zero=0.0):
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
    zero : float
        Fraction of zero-length intervals.

    Returns
    -------
    polars.DataFrame
        String columns chrom, start, end, name, fdr.
    """
    start = rng.integers(0, span, n)
    length = rng.integers(1, max_len, n)
    length[rng.random(n) < zero] = 0
    return pl.DataFrame(
        {
            "chrom": rng.choice(CHROMS, n).tolist(),
            "start": start.astype(str).tolist(),
            "end": (start + length).astype(str).tolist(),
            "name": [f"r{i}" for i in range(n)],
            "fdr": rng.choice(FDRS, n).tolist(),
        }
    )


def arrays(df):
    """Return (codes, starts, ends) of a table with one fixed chrom numbering.

    Parameters
    ----------
    df : polars.DataFrame
        Table with chrom, start and end columns.

    Returns
    -------
    tuple of numpy.ndarray
        Chrom codes, starts and ends.
    """
    codes = np.array([CHROMS.index(c) for c in df["chrom"]], dtype=np.int64)
    return (
        codes,
        df["start"].cast(pl.Int64).to_numpy(),
        df["end"].cast(pl.Int64).to_numpy(),
    )


def write(df, tmp_path, name):
    """Write a table and return its path.

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
    _bed.write_tsv(df, path)
    return path


def text(df):
    """Render a table as file text.

    Parameters
    ----------
    df : polars.DataFrame
        String columns.

    Returns
    -------
    str
        Tab-separated lines, each ending in a newline.
    """
    return "".join(f"{line}\n" for line in _bed.lines(df).to_list())


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("distance", [0, 1])
def test_merge(bedtools, tmp_path, seed, distance):
    """merge_clusters reproduces `bedtools merge -d` with collapse and min."""
    rng = np.random.default_rng(seed)
    df = random_bed(rng, 2000, span=20000, max_len=200)
    codes, starts, ends = arrays(df)
    order = np.lexsort((starts, codes))
    df, codes, starts, ends = df[order], codes[order], starts[order], ends[order]
    path = write(df, tmp_path, "a.bed")

    cluster = _bed.merge_clusters(codes, starts, ends, distance)
    ours = (
        df.with_columns(
            cluster=cluster,
            end=pl.Series(ends),
            fdr=pl.Series([float(x) for x in df["fdr"]]),
        )
        .group_by("cluster", maintain_order=True)
        .agg(
            pl.col("chrom").first(),
            pl.col("start").first(),
            pl.col("end").max().cast(pl.String),
            pl.col("name").str.join(","),
            pl.col("fdr").min(),
        )
        .drop("cluster")
    )
    ours = ours.with_columns(fdr=pl.Series([f"{v:.10g}" for v in ours["fdr"]]))
    expected = bedtools(
        "merge", "-i", path, "-d", distance, "-c", "4,5", "-o", "collapse,min"
    )
    assert text(ours) == expected


@pytest.mark.parametrize("seed", range(4))
def test_overlaps_any(bedtools, tmp_path, seed):
    """overlaps_any reproduces `bedtools intersect -u` and `-v`, zero-length included."""
    rng = np.random.default_rng(seed)
    a = random_bed(rng, 1500, zero=0.1)
    b = random_bed(rng, 1500, max_len=60, zero=0.2)
    pa, pb = write(a, tmp_path, "a.bed"), write(b, tmp_path, "b.bed")
    hit = _bed.overlaps_any(arrays(a), arrays(b))
    assert text(a.filter(hit)) == bedtools("intersect", "-u", "-a", pa, "-b", pb)
    assert text(a.filter(~hit)) == bedtools("intersect", "-v", "-a", pa, "-b", pb)


@pytest.mark.parametrize("seed", range(4))
def test_contains_any(bedtools, tmp_path, seed):
    """contains_any reproduces `bedtools intersect -F 1 -u` and `-v`."""
    rng = np.random.default_rng(seed)
    a = random_bed(rng, 1500, zero=0.1)
    b = random_bed(rng, 1500, max_len=60, zero=0.2)
    pa, pb = write(a, tmp_path, "a.bed"), write(b, tmp_path, "b.bed")
    hit = _bed.contains_any(arrays(a), arrays(b))
    assert text(a.filter(hit)) == bedtools(
        "intersect", "-F", 1, "-u", "-a", pa, "-b", pb
    )
    assert text(a.filter(~hit)) == bedtools(
        "intersect", "-F", 1, "-v", "-a", pa, "-b", pb
    )


@pytest.mark.parametrize("seed", range(4))
def test_overlap_pairs(bedtools, tmp_path, seed):
    """overlap_pairs gives the pairs and overlaps of `bedtools intersect -wo`."""
    rng = np.random.default_rng(seed)
    a = random_bed(rng, 1000)
    b = random_bed(rng, 1000, max_len=60, zero=0.2)
    pa, pb = write(a, tmp_path, "a.bed"), write(b, tmp_path, "b.bed")
    ia, ib, overlap = _bed.overlap_pairs(arrays(a), arrays(b))
    a_lines, b_lines = _bed.lines(a).to_list(), _bed.lines(b).to_list()
    ours = [f"{a_lines[i]}\t{b_lines[j]}\t{o}" for i, j, o in zip(ia, ib, overlap)]
    expected = bedtools("intersect", "-wo", "-a", pa, "-b", pb).splitlines()
    assert sorted(ours) == sorted(expected)
    # Hits come grouped by A record, in A order.
    assert [line.split("\t")[3] for line in ours] == [
        line.split("\t")[3] for line in expected
    ]


def test_overlap_pairs_rejects_zero_length_a():
    """Zero-length A intervals raise instead of mimicking bedtools' print quirk."""
    a = (np.array([0]), np.array([5]), np.array([5]))
    with pytest.raises(ValueError):
        _bed.overlap_pairs(a, a)


def tied_lines(rng, n):
    """Draw BED-like lines with many (chrom, start) ties and mixed-case names.

    Parameters
    ----------
    rng : numpy.random.Generator
        Random source.
    n : int
        Number of lines.

    Returns
    -------
    polars.DataFrame
        String columns chrom, start, end, name, value.
    """
    start = rng.integers(0, 40, n)
    names = [
        f"ENCFF{rng.integers(0, 5):03d}ABC-{rng.integers(1, 2000)}" for _ in range(n)
    ]
    return pl.DataFrame(
        {
            "chrom": rng.choice(CHROMS, n).tolist(),
            "start": start.astype(str).tolist(),
            "end": (start + rng.integers(0, 2000, n)).astype(str).tolist(),
            "name": [s.lower() if rng.random() < 0.2 else s for s in names],
            "value": rng.choice(FDRS, n).tolist(),
        }
    )


@pytest.mark.parametrize("seed", range(3))
def test_gnu_sort_bed(gnu_sort, seed):
    """gnu_sort reproduces `sort -k1,1 -k2,2n` with its whole-line tie-break."""
    df = tied_lines(np.random.default_rng(seed), 3000)
    order = _bed.gnu_sort(df, text_keys=["chrom"], numeric_keys=["start"])
    assert text(df[order]) == gnu_sort(["-k1,1", "-k2,2n"], text(df))


@pytest.mark.parametrize("seed", range(3))
def test_gnu_sort_field_key(gnu_sort, seed):
    """A `-k4,4` key is field 4 with its leading tab."""
    df = tied_lines(np.random.default_rng(seed), 3000)
    keyed = df.with_columns(key="\t" + pl.col("name"))
    order = _bed.gnu_sort(keyed, text_keys=["key"])
    assert text(df[order]) == gnu_sort(["-k4,4"], text(df))


@pytest.mark.parametrize("seed", range(3))
def test_sort_unique(gnu_sort, seed):
    """sort_unique reproduces `sort -u`."""
    df = tied_lines(np.random.default_rng(seed), 3000)
    line_list = _bed.lines(df).to_list()
    line_list += line_list[::7]
    expected = gnu_sort(["-u"], "".join(f"{s}\n" for s in line_list))
    assert "".join(f"{s}\n" for s in _bed.sort_unique(line_list)) == expected


def test_read_tsv_roundtrip(tmp_path):
    """read_tsv and write_tsv keep every byte, gzip input included."""
    content = "chr1\t10\t20\tn\t1.0e-05\t\nchr2\t5\t7\tm\t0.0468394\tx\n"
    plain = tmp_path / "a.tsv"
    plain.write_text(content)
    df = _bed.read_tsv(plain, ["c", "s", "e", "n", "f", "x"])
    out = tmp_path / "b.tsv"
    _bed.write_tsv(df, out)
    assert out.read_text() == content
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    assert _bed.read_tsv(empty, ["c", "s"]).shape == (0, 2)


@pytest.mark.parametrize("seed", range(4))
def test_closest_pairs(bedtools, tmp_path, seed):
    """closest_pairs lists the pairs of `bedtools closest -d`, ties and order included."""
    rng = np.random.default_rng(seed)
    a = random_bed(rng, 800, span=20000, max_len=300)
    b = random_bed(rng, 400, span=20000, max_len=40, zero=0.5)
    # Repeat some B rows so that several B share a position.
    b = pl.concat([b, b.sample(150, seed=seed)])
    a_codes, a_starts, _ = arrays(a)
    b_codes, b_starts, _ = arrays(b)
    a = a[np.lexsort((a_starts, a_codes))]
    b = b[np.lexsort((b_starts, b_codes))]
    pa, pb = write(a, tmp_path, "a.bed"), write(b, tmp_path, "b.bed")
    # Number B's chromosomes against its row order, as shared_codes may.
    b_arrays = arrays(b)
    b_arrays = (len(CHROMS) - 1 - b_arrays[0], *b_arrays[1:])
    a_arrays = arrays(a)
    a_arrays = (len(CHROMS) - 1 - a_arrays[0], *a_arrays[1:])
    ia, ib = _bed.closest_pairs(a_arrays, b_arrays)
    a_lines, b_lines = _bed.lines(a).to_list(), _bed.lines(b).to_list()
    ours = [f"{a_lines[i]}\t{b_lines[j]}" for i, j in zip(ia, ib)]
    expected = [
        line.rsplit("\t", 1)[0]
        for line in bedtools("closest", "-d", "-a", pa, "-b", pb).splitlines()
    ]
    assert ours == expected
