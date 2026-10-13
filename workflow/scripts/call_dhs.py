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
    """Return the order filter-long-double.py reads from an FDR.

    Parameters
    ----------
    text : str
        FDR as written in the track.

    Returns
    -------
    int
        Exponent after the first '-' of the text, printed as '%.2E' first if it
        has no exponent. `KEEP` if the printed text has a '+'.
    """
    if "e" not in text and "E" not in text:
        text = "%.2E" % Decimal(text)  # noqa: UP031 (the original's formatting)
    if "+" in text:
        return KEEP
    return int(text.split("-")[1])


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

    def merged(self, rows):
        """Merge positions like `bedtools merge -d 1 -c 5 -o min`.

        Parameters
        ----------
        rows : numpy.ndarray
            Indices of the kept positions, increasing.

        Returns
        -------
        tuple of (dict of str to numpy.ndarray, numpy.ndarray)
            Code, start, end and min FDR of each merged interval, in order, and
            the interval index of each row.
        """
        codes, starts, ends = self.codes[rows], self.starts[rows], self.ends[rows]
        cluster = _bed.merge_clusters(codes, starts, ends, GAP)
        first = np.flatnonzero(np.diff(cluster, prepend=-1))
        merged = {
            "code": codes[first],
            "start": starts[first],
            "end": np.maximum.reduceat(ends, first) if len(rows) else ends,
            "fdr": np.minimum.reduceat(self.fdrs[rows], first) if len(rows) else ends,
        }
        return merged, cluster


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
    tuple of (list of dict, dict)
        DHSs in calling order, threshold by threshold, and the excluded
        intervals.
    """
    rows = np.arange(len(track.orders))
    changes = np.unique(track.orders[(track.orders >= 3) & (track.orders <= min_p)])
    dhs = []
    for j in [2, *changes.tolist()]:
        rows = rows[track.orders[rows] > j]
        merged, cluster = track.merged(rows)
        length = merged["end"] - merged["start"]
        dhs.append(_select(merged, (length >= MIN_LENGTH) & (length <= MAX_DHS_LENGTH)))
        wide = length > MAX_DHS_LENGTH
        rows = rows[wide[cluster]]
        excluded = _select(merged, wide)
        if not len(rows):
            break
    return dhs, excluded


def _select(intervals, mask):
    """Keep the intervals where `mask` is true.

    Parameters
    ----------
    intervals : dict of str to numpy.ndarray
        Interval arrays.
    mask : numpy.ndarray
        Flag per interval.

    Returns
    -------
    dict of str to numpy.ndarray
        Selected intervals.
    """
    return {k: v[mask] for k, v in intervals.items()}


def write(track, parts, path):
    """Write intervals as chrom, start, end and min FDR, printed as bedtools does.

    Parameters
    ----------
    track : Track
        Source of the chromosome names.
    parts : list of dict
        Intervals, written in order.
    path : str
        BED to write.
    """
    with open(path, "w") as f:
        for part in parts:
            rows = zip(
                part["code"].tolist(),
                part["start"].tolist(),
                part["end"].tolist(),
                part["fdr"].tolist(),
            )
            # bedtools prints with std::setprecision(10).
            f.writelines(
                f"{track.chroms[code]}\t{start}\t{end}\t{fdr:.10g}\n"
                for code, start, end, fdr in rows
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
    write(track, [excluded], args.excluded)


if __name__ == "__main__":
    main()
