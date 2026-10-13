"""Compare each pipeline output with its May 2022 reference file.

Select a step with `-m <step>` and a genome with `-k <genome>`.
"""

import gzip
from pathlib import Path

# Columns 5-7 hold ;-joined sets. The original printed them in set() order.
LIST_COLUMNS = (4, 5, 6)
# Bound on z-score gaps from glibc's FMA versus non-FMA log. Observed: 9e-16.
Z_TOLERANCE = 1e-14


def read_ours(path):
    """Open our file as a byte stream, decompressing gzip on the fly.

    Parameters
    ----------
    path : pathlib.Path
        Our file. A `.gz` suffix marks gzip or bgzip.

    Returns
    -------
    io.BufferedIOBase
        Stream of the uncompressed bytes.
    """
    return gzip.open(path, "rb") if path.suffix == ".gz" else open(path, "rb")


def same_bytes(a, b):
    """Return whether our file holds the bytes of the reference file.

    Parameters
    ----------
    a, b : pathlib.Path
        Our file, read through `read_ours`, and the reference file.

    Returns
    -------
    bool
        True if the contents are equal.
    """
    with read_ours(a) as fa, open(b, "rb") as fb:
        while (chunk := fa.read(1 << 24)) == fb.read(1 << 24):
            if not chunk:
                return True
    return False


def first_difference(a, b, normalize_a, normalize_b):
    """Find the first line where two text files differ.

    Parameters
    ----------
    a, b : pathlib.Path
        Our file, read through `read_ours`, and the reference file.
    normalize_a, normalize_b : callable
        Map a line of `a` or `b` to its comparable form.

    Returns
    -------
    str or None
        Description of the first difference, None if there is none.
    """
    with read_ours(a) as fa, open(b, "rb") as fb:
        n = 0
        for n, (la, lb) in enumerate(zip(fa, fb), start=1):
            if normalize_a(la) != normalize_b(lb):
                return f"line {n}: ours {la[:300]!r} reference {lb[:300]!r}"
        rest_a, rest_b = fa.readline(), fb.readline()
    if rest_a or rest_b:
        return f"{'ours' if rest_a else 'reference'} has more lines after line {n}"
    return None


def strip(line):
    """Drop the newline of a line.

    Parameters
    ----------
    line : bytes
        One line.

    Returns
    -------
    bytes
        The line without its trailing newline.
    """
    return line.rstrip(b"\n")


def first_8_columns(line):
    """Keep the first 8 tab-separated fields of a line.

    Parameters
    ----------
    line : bytes
        One line.

    Returns
    -------
    bytes
        Fields 1-8, tab-joined, without the newline.
    """
    return b"\t".join(strip(line).split(b"\t")[:8])


def except_column_6(line):
    """Drop field 6 of a line.

    Parameters
    ----------
    line : bytes
        One line.

    Returns
    -------
    list of bytes
        Fields other than field 6.
    """
    fields = strip(line).split(b"\t")
    return fields[:5] + fields[6:]


def z_close(ours, ref):
    """Compare z-score tables, allowing rounding-level gaps in column 2.

    Parameters
    ----------
    ours, ref : pathlib.Path
        Signal or maxZ tables, ours first.

    Returns
    -------
    str or None
        First difference beyond the tolerance, None if there is none.
    """
    with open(ours, "rb") as fa, open(ref, "rb") as fb:
        n = 0
        for n, (la, lb) in enumerate(zip(fa, fb), start=1):
            a, b = strip(la).split(b"\t"), strip(lb).split(b"\t")
            za, zb = float(a[1]), float(b[1])
            if a[:1] + a[2:] != b[:1] + b[2:] or abs(za - zb) > Z_TOLERANCE:
                return f"line {n}: ours {la[:300]!r} reference {lb[:300]!r}"
        rest_a, rest_b = fa.readline(), fb.readline()
    if rest_a or rest_b:
        return f"{'ours' if rest_a else 'reference'} has more lines after line {n}"
    return None


def sort_lists(line):
    """Sort the ;-joined items of the set-valued columns of a line.

    Parameters
    ----------
    line : bytes
        One line.

    Returns
    -------
    list of bytes
        Fields, with the set-valued ones in sorted item order.
    """
    fields = strip(line).split(b"\t")
    for i in LIST_COLUMNS:
        if i < len(fields):
            fields[i] = b";".join(sorted(fields[i].split(b";")))
    return fields


def difference(ours, ref, mode):
    """Describe how our file differs from the reference.

    Parameters
    ----------
    ours, ref : pathlib.Path
        Our file and the reference file.
    mode : str
        `bytes`, `ref-cols-1-8`, `except-col-6`, `z-close` or `sorted-lists`.

    Returns
    -------
    str or None
        First difference, None if the files match under `mode`.
    """
    if mode == "bytes":
        if same_bytes(ours, ref):
            return None
        return first_difference(ours, ref, strip, strip) or "same lines, other bytes"
    if mode == "ref-cols-1-8":
        return first_difference(ours, ref, strip, first_8_columns)
    if mode == "z-close":
        return z_close(ours, ref)
    if mode == "except-col-6":
        return first_difference(ours, ref, except_column_6, except_column_6)
    return first_difference(ours, ref, sort_lists, sort_lists)


def test_matches_reference(case):
    """Our output equals its reference file under the case's mode.

    Parameters
    ----------
    case : tuple of (pathlib.Path, pathlib.Path, str)
        Our file, the reference file and the comparison mode.
    """
    ours, ref, mode = case
    assert Path(ref).exists(), f"no reference file {ref}"
    detail = difference(ours, ref, mode)
    assert detail is None, detail
