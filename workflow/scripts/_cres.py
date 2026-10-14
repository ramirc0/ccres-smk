"""Pieces shared by the cell type-agnostic (step 6) and specific (step 7) calls."""

import _bed
import numpy as np
import polars as pl

HIGH = 1.64
CENTER = 200
FAR = 2000
TSS = ["chrom", "start", "end", "transcript", "score", "strand", "gene"]


def select(keys, items, negate=False, field=0):
    """Keep the items whose key is in `keys`, or with `negate` not in it.

    Port of `awk 'FNR==NR {x[$k];next} (key in x)' first second`. awk read an
    empty first file as no file and printed nothing, even under `!(key in x)`.
    Here no keys select every item under `negate`.

    Parameters
    ----------
    keys : list of str
        Keys of the first file.
    items : list
        Items of the second file: all keys, or all sequences such as tuples
        and rows that hold the key at `field`.
    negate : bool
        Select the items whose key is not in `keys`.
    field : int
        Position of the key in a sequence item.

    Returns
    -------
    list
        Selected items, in order.
    """
    keys = set(keys)
    if items and isinstance(items[0], str):
        return [s for s in items if (s in keys) != negate]
    return [s for s in items if (s[field] in keys) != negate]


def above(lines):
    """Keep the signal lines whose value exceeds 1.64.

    Parameters
    ----------
    lines : list of tuple of (str, str)
        Name and value text of each line.

    Returns
    -------
    list of tuple of (str, str)
        Lines above 1.64.
    """
    return [line for line in lines if float(line[1]) > HIGH]


def at_most(lines):
    """Keep the signal lines whose value is at most 1.64.

    Parameters
    ----------
    lines : list of tuple of (str, str)
        Name and value text of each line.

    Returns
    -------
    list of tuple of (str, str)
        Lines at or below 1.64.
    """
    return [line for line in lines if float(line[1]) <= HIGH]


def read_signal(path):
    """Read the name and value columns of a signal or maxZ table.

    Parameters
    ----------
    path : str or None
        Table path, or None for a mark the biosample lacks.

    Returns
    -------
    list of tuple of (str, str)
        Name and value text of each line, none without a table.
    """
    if path is None:
        return []
    with open(path) as f:
        return [tuple(line.split("\t")[:2]) for line in f]


def tss_groups(names, tss, tss_iv, window_iv, region_iv):
    """Split regions into TSS, proximal and distal groups.

    Port of the shared block of 6_Classify-cCREs.sh and Split-cCREs.*.sh,
    with calculate-center-distance.py. TSS regions overlap a TSS. Proximal
    regions overlap a TSS +-2 kb window. Distal regions overlap none. A
    proximal region whose center lies within 200 bp of its first closest TSS
    moves to TSS. One over 2 kb away moves to distal.

    Parameters
    ----------
    names : list of str
        Name of each region, in file order. The key of the selections.
    tss : polars.DataFrame
        TSS table.
    tss_iv, window_iv, region_iv : tuple of numpy.ndarray
        Intervals of the TSSs, their windows and the regions, on one
        chromosome numbering.

    Returns
    -------
    dict of str to list of str
        Region names of the tss, prox and distal groups, in file order.
    """
    on_tss = _bed.overlaps_any(region_iv, tss_iv)
    in_window = _bed.overlaps_any(region_iv, window_iv)
    names = np.array(names, dtype=object)
    tss_names = names[on_tss].tolist()
    prox_rows = np.flatnonzero(~on_tss & in_window)
    prox_names = names[prox_rows].tolist()
    distal_names = names[~in_window].tolist()

    prox_iv = tuple(a[prox_rows] for a in region_iv)
    ia, ib = _bed.closest_pairs(prox_iv, tss_iv)
    first = np.flatnonzero(np.diff(ia, prepend=-1))
    site = ib[first]
    center = (prox_iv[1][ia[first]] + prox_iv[2][ia[first]]) / 2
    start = tss["start"].cast(pl.Int64).to_numpy()[site]
    plus = tss["strand"].to_numpy()[site] == "+"
    distance = np.where(plus, start - center, center - start)
    keyed = [prox_names[i] for i in ia[first]]
    near = [n for n, d in zip(keyed, distance) if -CENTER <= d <= CENTER]
    far = [n for n, d in zip(keyed, distance) if d < -FAR or d > FAR]

    return {
        "tss": tss_names + select(near, prox_names),
        "prox": select(near + far, prox_names, negate=True),
        "distal": distal_names + select(far, prox_names),
    }
