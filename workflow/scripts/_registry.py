"""Representative-peak selection (steps 2 and 3) and region accessioning."""

import _bed
import numpy as np

GENOME_ID = {"hg19": "EH37", "hg38": "EH38", "mm10": "EM10"}
MODE_ID = {"rDHS": "D", "cCRE": "E", "TF": "F"}


def representative_peaks(codes, starts, ends, signal, names):
    """Pick representative peaks like the merge, pick and intersect loop.

    Each pass merges the remaining peaks (`bedtools merge`, book-ended
    included), picks the first peak of highest signal per cluster
    (`pick-best-peak.py`), keeps every remaining peak with a picked name
    (the awk join) and drops every remaining peak that overlaps a pick
    (`bedtools intersect -v`). Passes repeat until no peak remains.

    Parameters
    ----------
    codes, starts, ends : numpy.ndarray
        Chromosome codes and coordinates, in `sort -k1,1 -k2,2n` order.
    signal : numpy.ndarray
        Signal as float64, parsed as Python's float() parses the text.
    names : numpy.ndarray
        Integer code of each peak's name. Equal names share a code.

    Returns
    -------
    numpy.ndarray
        Row indices of the representative peaks, pass by pass, in row order
        within a pass.
    """
    if np.isnan(signal).any():
        raise ValueError("NaN signal: Python's max() would depend on its position.")
    unique_names = len(np.unique(names)) == len(names)
    remaining = np.arange(len(codes))
    picks = []
    while len(remaining):
        c, s, e = codes[remaining], starts[remaining], ends[remaining]
        cluster = _bed.merge_clusters(c, s, e)
        first = np.flatnonzero(np.diff(cluster, prepend=-1))
        sig = signal[remaining]
        top = np.maximum.reduceat(sig, first)[cluster]
        position = np.where(sig == top, np.arange(len(sig)), len(sig))
        best = np.minimum.reduceat(position, first)
        if unique_names:
            chosen = best
        else:
            chosen = np.flatnonzero(np.isin(names[remaining], names[remaining[best]]))
        picks.append(remaining[chosen])
        hit = _bed.overlaps_any((c, s, e), (c[chosen], s[chosen], e[chosen]))
        remaining = remaining[~hit]
    return np.concatenate(picks) if picks else np.array([], dtype=np.int64)


def previous_accessions(rows, mode):
    """Index the regions of a previous release, like make-region-accession.py.

    Parameters
    ----------
    rows : iterable of sequence of str
        Fields of each previous region: chrom, start, end, accession.
    mode : str
        `rDHS`, `cCRE` or `TF`.

    Returns
    -------
    tuple of (dict of str to str, int)
        Accession per concatenated `chrom + start + end` key, last row winning,
        and the largest accession number.
    """
    previous, top = {}, 0
    for chrom, start, end, accession, *_ in rows:
        previous[chrom + start + end] = accession
        top = max(top, int(accession.split(MODE_ID[mode])[-1]))
    return previous, top


def accessions(regions, genome, mode, previous=None):
    """Assign accessions to regions like make-region-accession.py.

    A region whose `chrom + start + end` string matches a previous region
    reuses that accession, with EH37 renamed EH38. The key has no separators,
    so distinct coordinates can collide. That quirk is kept. Other regions
    get new numbers above the previous maximum, in input order.

    Parameters
    ----------
    regions : iterable of tuple of str
        `(chrom, start, end)` of each region, in output order.
    genome : str
        `hg19`, `hg38` or `mm10`.
    mode : str
        `rDHS`, `cCRE` or `TF`.
    previous : tuple of (dict of str to str, int), optional
        Output of `previous_accessions`.

    Returns
    -------
    list of str
        Accession of each region.
    """
    known, top = previous or ({}, 0)
    prefix = GENOME_ID[genome] + MODE_ID[mode]
    out = []
    i = top + 1
    for chrom, start, end in regions:
        key = chrom + start + end
        if key in known:
            out.append(known[key].replace("EH37", "EH38"))
        else:
            out.append(f"{prefix}{i:07d}")
            i += 1
    return out
