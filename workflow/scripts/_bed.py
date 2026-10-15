"""Exact BED text IO, GNU sort under the C locale, and bedtools 2.30.0 interval ops.

Tables are polars DataFrames of string columns, read and written byte for byte.
Interval ops take numpy arrays and reproduce bedtools 2.30.0. bedtools reads a
zero-length interval [s, s) as [s - 1, s + 1), and so do these ops.

Text sorts by bytes, as `LC_ALL=C sort` does. The release sorted under
en_US.UTF-8, which ignores blanks and punctuation on its first pass.
"""

import numpy as np
import polars as pl

# Coordinates of one chromosome fit in 2**34. Adding chrom * 2**34 makes one
# sorted axis of (chrom, position) pairs.
_CHROM_SHIFT = 34


def read_tsv(path, columns):
    """Read a tab-separated file as string columns, exactly as written.

    Parameters
    ----------
    path : str
        File to read. Gzip input is decompressed.
    columns : list of str
        Column names. Their number must match the file.

    Returns
    -------
    polars.DataFrame
        One string column per field. An empty file gives zero rows.
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
        return pl.DataFrame(schema=dict.fromkeys(columns, pl.String))
    if df.width != len(columns):
        raise ValueError(f"{path}: {df.width} columns, expected {len(columns)}.")
    # read_csv returns thousands of chunks per column, which make gathers slow.
    return df.rename(dict(zip(df.columns, columns))).rechunk()


def lines(df):
    """Join the columns of a table into tab-separated lines.

    Parameters
    ----------
    df : polars.DataFrame
        String columns.

    Returns
    -------
    polars.Series
        One line per row, without the newline.
    """
    return df.select(pl.concat_str(pl.all(), separator="\t").alias("line"))["line"]


def write_tsv(df, path):
    """Write a table as tab-separated lines.

    Parameters
    ----------
    df : polars.DataFrame
        String columns, written in order.
    path : str
        File to write.
    """
    lines(df).to_frame().write_csv(path, include_header=False, quote_style="never")


def gnu_sort(df, text_keys=(), numeric_keys=()):
    """Order rows like `LC_ALL=C sort` with these keys, without -u.

    Rows compare by the text keys, then the numeric keys, then by whole line,
    GNU sort's last resort. Text compares by bytes. No field holds a tab, so
    whole lines compare as their columns do, in order.

    Parameters
    ----------
    df : polars.DataFrame
        String columns. A row's line is its columns joined by tabs.
    text_keys : sequence of str
        Columns compared as text.
    numeric_keys : sequence of str
        Integer columns compared by value (`-n`).

    Returns
    -------
    numpy.ndarray
        Row indices in sorted order.
    """
    numbers = [pl.col(name).cast(pl.Int64, strict=True) for name in numeric_keys]
    return (
        df.with_row_index("_row")
        .sort([*text_keys, *numbers, *df.columns], maintain_order=True)
        .get_column("_row")
        .to_numpy()
    )


def sort_unique(line_list):
    """Sort lines like `LC_ALL=C sort -u`.

    Parameters
    ----------
    line_list : list of str
        Lines without newlines.

    Returns
    -------
    list of str
        Distinct lines in byte order.
    """
    return sorted(set(line_list))


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


def chrom_codes(chroms):
    """Code chromosomes by order of first appearance.

    Parameters
    ----------
    chroms : polars.Series
        Chromosome of each row.

    Returns
    -------
    numpy.ndarray
        Code of each row, int64. Rows grouped by chromosome get increasing codes.
    """
    first_seen = chroms.unique(maintain_order=True).to_list()
    codes = {chrom: i for i, chrom in enumerate(first_seen)}
    return chroms.replace_strict(codes, return_dtype=pl.Int64).to_numpy()


def shared_codes(*chroms):
    """Code the chromosomes of several tables with one numbering.

    Parameters
    ----------
    *chroms : polars.Series
        Chromosome of each row, one series per table.

    Returns
    -------
    list of numpy.ndarray
        Codes per table as int64, in the order given.
    """
    first_seen = pl.concat(chroms).unique(maintain_order=True).to_list()
    codes = {chrom: i for i, chrom in enumerate(first_seen)}
    return [c.replace_strict(codes, return_dtype=pl.Int64).to_numpy() for c in chroms]


def python_floats(series):
    """Parse number texts as Python's float() does.

    polars rounds correctly, as float(), C's strtod and awk do. Texts that
    float() accepts but polars rejects, such as " 1" or "1_0", raise instead.

    Parameters
    ----------
    series : polars.Series
        Number texts.

    Returns
    -------
    numpy.ndarray
        Values as float64.
    """
    return series.cast(pl.Float64, strict=True).to_numpy()


def axis(codes, positions):
    """Place (chromosome code, position) pairs on one integer axis.

    Parameters
    ----------
    codes : numpy.ndarray
        Chromosome codes.
    positions : numpy.ndarray
        Positions, below 2**34.

    Returns
    -------
    numpy.ndarray
        Positions as int64, ordered by code first.
    """
    return (codes.astype(np.int64) << _CHROM_SHIFT) + positions


def bedtools_coords(starts, ends):
    """Apply bedtools' zero-length adjustment, [s, s) -> [s - 1, s + 1).

    Parameters
    ----------
    starts, ends : numpy.ndarray
        BED coordinates, int64.

    Returns
    -------
    tuple of numpy.ndarray
        Adjusted starts and ends.
    """
    zero = starts == ends
    return starts - zero, ends + zero


def merge_clusters(codes, starts, ends, distance=0):
    """Assign rows to the clusters of `bedtools merge -d distance`.

    A row joins the open cluster if its start is at most the cluster's
    furthest end plus `distance`.

    Parameters
    ----------
    codes : numpy.ndarray
        Chromosome codes. Rows of a chromosome are contiguous, codes increase.
    starts, ends : numpy.ndarray
        Coordinates, int64, sorted by start within each chromosome.
        Zero-length intervals are rejected: bedtools widens them first.
    distance : int
        Largest gap that still merges.

    Returns
    -------
    numpy.ndarray
        Cluster index of each row, int64, from 0 in row order.
    """
    if np.any(starts == ends):
        raise ValueError("merge_clusters does not handle zero-length intervals.")
    s, e = axis(codes, starts), axis(codes, ends)
    if np.any(np.diff(s) < 0):
        raise ValueError("Rows are not sorted by chromosome and start.")
    new = np.ones(len(s), dtype=bool)
    new[1:] = s[1:] > np.maximum.accumulate(e)[:-1] + distance
    return np.cumsum(new) - 1


def _sorted_targets(b_codes, b_starts, b_ends):
    """Sort target intervals on the joint axis, with prefix maxima of ends.

    Parameters
    ----------
    b_codes, b_starts, b_ends : numpy.ndarray
        Target intervals, bedtools-adjusted.

    Returns
    -------
    tuple of numpy.ndarray
        Sorted axis starts, matching axis ends, prefix max of axis ends, and
        the original index of each sorted target.
    """
    bs, be = axis(b_codes, b_starts), axis(b_codes, b_ends)
    order = np.argsort(bs, kind="stable")
    bs, be = bs[order], be[order]
    return bs, be, np.maximum.accumulate(be) if len(be) else be, order


def overlaps_any(a, b):
    """Flag A intervals that overlap any B, as `bedtools intersect -u`.

    Parameters
    ----------
    a, b : tuple of numpy.ndarray
        `(codes, starts, ends)`. Codes MUST come from one shared numbering.

    Returns
    -------
    numpy.ndarray
        Flag per A interval.
    """
    a_s, a_e = bedtools_coords(a[1], a[2])
    b_s, b_e = bedtools_coords(b[1], b[2])
    bs, _, prefix_end, _ = _sorted_targets(b[0], b_s, b_e)
    as_, ae = axis(a[0], a_s), axis(a[0], a_e)
    # Targets starting before A ends. One of them overlaps iff its end is past A's start.
    n = np.searchsorted(bs, ae, side="left")
    hit = np.zeros(len(as_), dtype=bool)
    has = n > 0
    hit[has] = prefix_end[n[has] - 1] > as_[has]
    return hit


def contains_any(a, b):
    """Flag A intervals that cover all of some B, as `bedtools intersect -F 1 -u`.

    bedtools tests overlap / len(B) >= 1 in float32. Below 2**24 bp that is
    exact containment, which is all this function implements.

    Parameters
    ----------
    a, b : tuple of numpy.ndarray
        `(codes, starts, ends)`. Codes MUST come from one shared numbering.

    Returns
    -------
    numpy.ndarray
        Flag per A interval.
    """
    a_s, a_e = bedtools_coords(a[1], a[2])
    b_s, b_e = bedtools_coords(b[1], b[2])
    if len(b_e) and np.max(b_e - b_s) >= 1 << 24:
        raise ValueError("B intervals of 2**24 bp or more need float32 fractions.")
    bs, be, _, _ = _sorted_targets(b[0], b_s, b_e)
    as_, ae = axis(a[0], a_s), axis(a[0], a_e)
    # Targets starting inside A: [lo, hi). One is contained iff the smallest end <= A's end.
    lo = np.searchsorted(bs, as_, side="left")
    hi = np.searchsorted(bs, ae, side="left")
    return _range_min(be, lo, hi) <= ae


def _range_min(values, lo, hi):
    """Minimum of `values[lo:hi]` for many ranges, with a sparse table.

    Parameters
    ----------
    values : numpy.ndarray
        Values, int64.
    lo, hi : numpy.ndarray
        Range bounds. Empty ranges give the int64 maximum.

    Returns
    -------
    numpy.ndarray
        Minimum per range, int64.
    """
    out = np.full(len(lo), np.iinfo(np.int64).max)
    size = hi - lo
    nonempty = size > 0
    if not np.any(nonempty):
        return out
    levels = [values]
    while (1 << len(levels)) <= size.max():
        prev, step = levels[-1], 1 << (len(levels) - 1)
        levels.append(np.minimum(prev[:-step], prev[step:]))
    lo, hi, size = lo[nonempty], hi[nonempty], size[nonempty]
    k = np.floor(np.log2(size)).astype(np.int64)
    result = np.empty(len(lo), dtype=np.int64)
    for level in np.unique(k):
        sel = k == level
        table = levels[level]
        result[sel] = np.minimum(table[lo[sel]], table[hi[sel] - (1 << level)])
    out[nonempty] = result
    return out


def overlap_pairs(a, b):
    """List overlapping (A, B) pairs and their overlap, as `bedtools intersect -wo`.

    Pairs come grouped by A in A order, then in B start order. bedtools emits a
    record's hits in bin order instead, so callers MUST not depend on the order
    within an A record.

    Parameters
    ----------
    a, b : tuple of numpy.ndarray
        `(codes, starts, ends)`. Codes MUST come from one shared numbering.
        Zero-length A intervals are rejected. bedtools prints their second and
        later hits from unwidened coordinates.

    Returns
    -------
    tuple of numpy.ndarray
        A index, B index and printed overlap of each pair. For a zero-length B,
        bedtools prints the widened overlap minus 2.
    """
    if np.any(a[1] == a[2]):
        raise ValueError("overlap_pairs does not handle zero-length A intervals.")
    a_s, a_e = a[1], a[2]
    b_s, b_e = bedtools_coords(b[1], b[2])
    bs, be, _, b_index = _sorted_targets(b[0], b_s, b_e)
    as_, ae = axis(a[0], a_s), axis(a[0], a_e)
    longest = int(np.max(be - bs)) if len(bs) else 0
    # Candidates start within `longest` before A's start, and before A's end.
    lo = np.searchsorted(bs, as_ - longest, side="left")
    hi = np.searchsorted(bs, ae, side="left")
    counts = np.maximum(hi - lo, 0)
    ia = np.repeat(np.arange(len(as_)), counts)
    offsets = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
    ib = np.repeat(lo, counts) + offsets
    overlap = np.minimum(ae[ia], be[ib]) - np.maximum(as_[ia], bs[ib])
    keep = overlap > 0
    ia, ib = ia[keep], b_index[ib[keep]]
    return ia, ib, overlap[keep] - 2 * (b[1][ib] == b[2][ib])


def closest_pairs(a, b):
    """List the closest B of each A in `bedtools closest -t all` output order.

    Overlapping B intervals come first, at distance 0. Without one, the
    nearest B ending at or before A's start (left) and the nearest B starting
    at or after A's end (right) compete. The nearer side wins, and a tie
    reports left before right. Within a side, tied B intervals come in B row
    order. Coordinates are bedtools-adjusted.

    Parameters
    ----------
    a, b : tuple of numpy.ndarray
        `(codes, starts, ends)`. Codes MUST come from one shared numbering. B
        rows MUST be sorted by start within each chromosome, as bedtools
        requires. The chromosome order is free.

    Returns
    -------
    tuple of numpy.ndarray
        A index and B index of each reported pair, grouped by A in A order.
    """
    a_s, a_e = bedtools_coords(a[1], a[2])
    b_s, b_e = bedtools_coords(b[1], b[2])
    same_chrom = b[0][1:] == b[0][:-1]
    if np.any(b_s[1:][same_chrom] < b_s[:-1][same_chrom]):
        raise ValueError("B rows are not sorted by start within a chromosome.")
    as_, ae = axis(a[0], a_s), axis(a[0], a_e)
    bs, be = axis(b[0], b_s), axis(b[0], b_e)
    by_start = np.argsort(bs, kind="stable")
    starts = bs[by_start]
    n = len(as_)

    # Overlaps: B starting before A's end and ending after A's start.
    ov_a, ov_b, _ = overlap_pairs((a[0], a_s, a_e), (b[0], b_s, b_e))
    has_overlap = np.zeros(n, dtype=bool)
    has_overlap[ov_a] = True

    # Left: largest end <= A's start, on A's chromosome.
    by_end = np.lexsort((np.arange(len(be)), be))
    ends = be[by_end]
    left_end = np.full(n, -1, dtype=np.int64)
    k = np.searchsorted(ends, as_, side="right") - 1
    ok = (k >= 0) & ~has_overlap
    ok[ok] = (ends[k[ok]] >> _CHROM_SHIFT) == (as_[ok] >> _CHROM_SHIFT)
    left_end[ok] = ends[k[ok]]

    # Right: smallest start >= A's end, on A's chromosome.
    right_start = np.full(n, -1, dtype=np.int64)
    k = np.searchsorted(starts, ae, side="left")
    ok = (k < len(starts)) & ~has_overlap
    ok[ok] = (starts[k[ok]] >> _CHROM_SHIFT) == (ae[ok] >> _CHROM_SHIFT)
    right_start[ok] = starts[k[ok]]

    left_dist = np.where(left_end >= 0, as_ - left_end + 1, np.iinfo(np.int64).max)
    right_dist = np.where(
        right_start >= 0, right_start - ae + 1, np.iinfo(np.int64).max
    )
    use_left = (left_end >= 0) & (left_dist <= right_dist)
    use_right = (right_start >= 0) & (right_dist <= left_dist)

    pairs = []
    for i in range(n):
        if has_overlap[i]:
            continue
        if use_left[i]:
            lo = np.searchsorted(ends, left_end[i], side="left")
            hi = np.searchsorted(ends, left_end[i], side="right")
            pairs.extend((i, j) for j in by_end[lo:hi])
        if use_right[i]:
            lo = np.searchsorted(starts, right_start[i], side="left")
            hi = np.searchsorted(starts, right_start[i], side="right")
            pairs.extend((i, j) for j in by_start[lo:hi])
    near = np.array(pairs, dtype=np.int64).reshape(-1, 2)
    ia = np.concatenate([ov_a, near[:, 0]])
    ib = np.concatenate([ov_b, near[:, 1]])
    order = np.lexsort((np.arange(len(ia)), ia))
    return ia[order], ib[order]
