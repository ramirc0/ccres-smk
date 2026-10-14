"""Step 0: call DHSs from one Hotspot2 FDR track by the threshold sweep.

Port of 0_Call-DHSs.sh with filter-long-double.py. At each threshold 1E-j,
j = 2..min_p, the original kept the positions whose FDR order exceeds j, merged
them with `bedtools merge -d 1 -c 5 -o min` and dropped merged intervals under
50 bp. The 1E-2 intervals of at most 350 bp are DHSs. At each later threshold,
intervals of at most 350 bp that overlap no DHS yet join them. The 1E-min_p
intervals that overlap no DHS are the excluded ones.

A threshold changes the kept positions only where some position has exactly
that FDR order, so this port merges only at those thresholds. Elsewhere the
merged set repeats, and a repeated set adds no DHS.
"""

import argparse
from decimal import Decimal

import _bed
import numpy as np
import polars as pl
from numba import njit

KEEP = np.iinfo(np.int64).max  # order of an FDR that passes every threshold
GAP = 1  # bedtools merge -d
MIN_LENGTH = 50
MAX_DHS_LENGTH = 350


def build_parser():
    """Build the call_dhs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the call_dhs.py flags.
    """
    p = argparse.ArgumentParser(description="Call DHSs from a Hotspot2 FDR track.")
    p.add_argument("--hotspots", required=True, help="Hotspot2 BED, gzipped.")
    p.add_argument(
        "--min-p", type=int, required=True, help="Last threshold, 1E-<min_p>."
    )
    p.add_argument("--dhs", required=True, help="DHS BED to write.")
    p.add_argument("--excluded", required=True, help="Excluded BED to write.")
    return p


def fdr_order(text):
    """Return the order of an FDR, the k with 10**-k <= FDR < 10**-(k - 1).

    An FDR passes threshold 1E-j when its order exceeds j. filter-long-double.py
    printed an FDR without an exponent as '%.2E' before reading its order, which
    rounds 0.000999745 into order 3, and an FDR of 1 into one that passes every
    threshold.

    Parameters
    ----------
    text : str
        FDR as written in the track.

    Returns
    -------
    int
        The order, read exactly from the decimal text. `KEEP` for 0.
    """
    value = Decimal(text)
    return KEEP if value == 0 else -value.adjusted()


def per_value(series, function, dtype):
    """Apply a function once per distinct value of a string column.

    Parameters
    ----------
    series : polars.Series
        Strings.
    function : callable
        Maps one string to one value.
    dtype : polars.DataType
        Output dtype.

    Returns
    -------
    numpy.ndarray
        `function` of each element.
    """
    distinct = series.unique().to_list()
    mapping = dict(zip(distinct, map(function, distinct)))
    return series.replace_strict(mapping, return_dtype=dtype).to_numpy()


class Track:
    """Hotspot positions as arrays, in file order.

    Parameters
    ----------
    df : polars.DataFrame
        String columns chrom, start, end, name, fdr.
    """

    def __init__(self, df):
        self.chroms = df["chrom"].unique(maintain_order=True).to_numpy()
        self.codes = _bed.chrom_codes(df["chrom"])
        self.starts = df["start"].cast(pl.Int64).to_numpy()
        self.ends = df["end"].cast(pl.Int64).to_numpy()
        self.orders = per_value(df["fdr"], fdr_order, pl.Int64)
        self.fdrs = _bed.python_floats(df["fdr"])


@njit(cache=True)
def _sweep(codes, starts, ends, fdrs, orders, thresholds):
    """Run the threshold sweep over positions sorted by chromosome and start.

    At each threshold the positions above it merge like `bedtools merge -d 1
    -c 5 -o min`. `rows` holds the remaining positions and is compacted in
    place, so it stays in file order.

    Parameters
    ----------
    codes, starts, ends : numpy.ndarray
        Chromosome codes and coordinates, int64.
    fdrs : numpy.ndarray
        FDR of each position, float64.
    orders : numpy.ndarray
        FDR order of each position, int64.
    thresholds : numpy.ndarray
        Thresholds in sweep order, int64.

    Returns
    -------
    tuple of numpy.ndarray
        First row, end and min FDR of each DHS in calling order, then of each
        interval still over MAX_DHS_LENGTH after the last threshold.
    """
    n = len(codes)
    rows = np.arange(n)
    first = np.empty(n, np.int64)
    end = np.empty(n, np.int64)
    fdr = np.empty(n, np.float64)
    begin = np.empty(n + 1, np.int64)
    dhs_first = np.empty(1024, np.int64)
    dhs_end = np.empty(1024, np.int64)
    dhs_fdr = np.empty(1024, np.float64)
    called = 0
    count = 0
    for j in thresholds:
        m = 0
        count = 0
        reach = 0
        for r in rows:
            if orders[r] <= j:
                continue
            if starts[r] == ends[r]:
                raise ValueError("Zero-length position: bedtools would widen it.")
            if m > 0:
                previous = rows[m - 1]
                if codes[r] < codes[previous] or (
                    codes[r] == codes[previous] and starts[r] < starts[previous]
                ):
                    raise ValueError("Positions are not sorted.")
            if m == 0 or codes[r] != codes[first[count - 1]] or starts[r] > reach + GAP:
                first[count] = r
                end[count] = ends[r]
                fdr[count] = fdrs[r]
                begin[count] = m
                count += 1
                reach = ends[r]
            else:
                end[count - 1] = max(end[count - 1], ends[r])
                reach = max(reach, ends[r])
                # np.minimum propagates NaN.
                if fdrs[r] < fdr[count - 1] or fdrs[r] != fdrs[r]:
                    fdr[count - 1] = fdrs[r]
            rows[m] = r
            m += 1
        begin[count] = m
        kept = 0
        for i in range(count):
            length = end[i] - starts[first[i]]
            if MIN_LENGTH <= length <= MAX_DHS_LENGTH:
                if called == len(dhs_first):
                    dhs_first = np.concatenate((dhs_first, np.empty_like(dhs_first)))
                    dhs_end = np.concatenate((dhs_end, np.empty_like(dhs_end)))
                    dhs_fdr = np.concatenate((dhs_fdr, np.empty_like(dhs_fdr)))
                dhs_first[called] = first[i]
                dhs_end[called] = end[i]
                dhs_fdr[called] = fdr[i]
                called += 1
            elif length > MAX_DHS_LENGTH:
                size = begin[i + 1] - begin[i]
                rows[kept : kept + size] = rows[begin[i] : begin[i + 1]]
                kept += size
        rows = rows[:kept]
        if not kept:
            break
    wide = np.array(
        [i for i in range(count) if end[i] - starts[first[i]] > MAX_DHS_LENGTH],
        dtype=np.int64,
    )
    return (
        dhs_first[:called],
        dhs_end[:called],
        dhs_fdr[:called],
        first[wide],
        end[wide],
        fdr[wide],
    )


def sweep(track, min_p):
    """Run the threshold sweep.

    Merged intervals only split as the threshold rises. An interval of 50-350 bp
    becomes a DHS, and its pieces could only overlap it later. An interval under
    50 bp only has shorter pieces, which the original dropped. Both leave the
    sweep. So no remaining interval overlaps a DHS, and each remaining one of
    50-350 bp is a new DHS. The intervals still over 350 bp after the last
    threshold are the excluded ones.

    Parameters
    ----------
    track : Track
        Hotspot positions.
    min_p : int
        Last threshold, 1E-<min_p>.

    Returns
    -------
    tuple of (tuple of numpy.ndarray, tuple of numpy.ndarray)
        First row, end and min FDR of the DHSs in calling order, threshold by
        threshold, and of the excluded intervals.
    """
    changes = np.unique(track.orders[(track.orders >= 3) & (track.orders <= min_p)])
    found = _sweep(
        track.codes,
        track.starts,
        track.ends,
        track.fdrs,
        track.orders,
        np.concatenate(([2], changes)),
    )
    return found[:3], found[3:]


def write(track, intervals, path):
    """Write intervals as chrom, start, end and min FDR, printed as bedtools does.

    Parameters
    ----------
    track : Track
        Source of the chromosome names and starts.
    intervals : tuple of numpy.ndarray
        First row, end and min FDR of each interval, written in order.
    path : str
        BED to write.
    """
    first, end, fdr = intervals
    rows = zip(
        track.codes[first].tolist(),
        track.starts[first].tolist(),
        end.tolist(),
        fdr.tolist(),
    )
    with open(path, "w") as f:
        # bedtools prints with std::setprecision(10).
        f.writelines(
            f"{track.chroms[code]}\t{start}\t{stop}\t{value:.10g}\n"
            for code, start, stop, value in rows
        )


def main(argv=None):
    """Call the DHSs and excluded intervals of one track.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    track = Track(
        _bed.read_tsv(args.hotspots, ["chrom", "start", "end", "name", "fdr"])
    )
    dhs, excluded = sweep(track, args.min_p)
    write(track, dhs, args.dhs)
    write(track, excluded, args.excluded)


if __name__ == "__main__":
    main()
