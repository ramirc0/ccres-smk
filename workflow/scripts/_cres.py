"""Pieces shared by the cell type-agnostic (step 6) and specific (step 7) calls."""

import _bed
import numpy as np
import polars as pl

HIGH = 1.64
CENTER = 200
FAR = 2000
TSS = ["chrom", "start", "end", "transcript", "score", "strand", "gene"]


def awk_in(first, second, negate=False, field=0):
    """Select like `awk 'FNR==NR {x[$k];next} (key in x)' first second`.

    awk reads an empty first file as no file, so it takes the second file as
    the first and prints nothing. That quirk is kept: it decides some classes
    when a signal file is empty.

    Parameters
    ----------
    first : list of str
        Keys of the first file, in file order.
    second : list
        Items of the second file: all keys, or all sequences such as tuples
        and rows that hold the key at `field`.
    negate : bool
        Select the items whose key is not in `first`, as `!(key in x)`.
    field : int
        Position of the key in a sequence item.

    Returns
    -------
    list
        Selected items of `second`, in order.
    """
    if not first:
        return []
    keys = set(first)
    if second and isinstance(second[0], str):
        return [s for s in second if (s in keys) != negate]
    return [s for s in second if (s[field] in keys) != negate]


def above(lines):
    """Keep signal lines whose value passes awk's `$2 > 1.64`.

    Parameters
    ----------
    lines : list of tuple of (str, str)
        Name and value text of each line.

    Returns
    -------
    list of tuple of (str, str)
        Lines above 1.64. A non-numeric value compares as text with "1.64",
        so the placeholder line ("blank", "blank") passes.
    """
    return [line for line in lines if _awk_compare(line[1]) > 0]


def at_most(lines):
    """Keep signal lines whose value passes awk's `$2 <= 1.64`.

    Parameters
    ----------
    lines : list of tuple of (str, str)
        Name and value text of each line.

    Returns
    -------
    list of tuple of (str, str)
        Lines at or below 1.64.
    """
    return [line for line in lines if _awk_compare(line[1]) <= 0]


def _awk_compare(text):
    """Compare a field with 1.64 as awk does.

    Parameters
    ----------
    text : str
        Field text.

    Returns
    -------
    int
        -1, 0 or 1. A field that does not look numeric compares as text with
        the constant printed by CONVFMT, "1.64".
    """
    try:
        a, b = float(text), HIGH
    except ValueError:
        a, b = text, f"{HIGH:.6g}"
    return (a > b) - (a < b)


def read_signal(path):
    """Read the name and value columns of a signal or maxZ table.

    Parameters
    ----------
    path : str or None
        Table path, or None for the placeholder blank.txt, whose single line
        is "blank blank".

    Returns
    -------
    list of tuple of (str, str)
        Name and value text of each line.
    """
    if path is None:
        return [("blank", "blank")]
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
        Name of each region, in file order. The key of the awk joins.
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
        "tss": tss_names + awk_in(near, prox_names),
        "prox": awk_in(near + far, prox_names, negate=True),
        "distal": distal_names + awk_in(far, prox_names),
    }
