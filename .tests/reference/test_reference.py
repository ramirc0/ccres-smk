"""Compare each pipeline output with its May 2022 reference file.

Select a step with `-m <step>` and a genome with `-k <genome>`.
"""

import subprocess
from pathlib import Path

# Columns 5-7 hold ;-joined sets. The original printed them in set() order.
LIST_COLUMNS = (4, 5, 6)


def same_bytes(a, b):
    """Return whether two files have identical bytes.

    Parameters
    ----------
    a, b : pathlib.Path
        Files to compare.

    Returns
    -------
    bool
        True if the contents are equal.
    """
    if a.stat().st_size != b.stat().st_size:
        return False
    with open(a, "rb") as fa, open(b, "rb") as fb:
        while chunk := fa.read(1 << 24):
            if chunk != fb.read(1 << 24):
                return False
    return True


def first_difference(a, b, normalize_a, normalize_b):
    """Find the first line where two text files differ.

    Parameters
    ----------
    a, b : pathlib.Path
        Files to compare, ours first.
    normalize_a, normalize_b : callable
        Map a line of `a` or `b` to its comparable form.

    Returns
    -------
    str or None
        Description of the first difference, None if there is none.
    """
    with open(a, "rb") as fa, open(b, "rb") as fb:
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


def bigbed_text(path, out):
    """Convert a bigBed to BED text with bigBedToBed.

    Parameters
    ----------
    path : pathlib.Path
        bigBed file.
    out : pathlib.Path
        BED text file to write.

    Returns
    -------
    pathlib.Path
        `out`.
    """
    subprocess.run(["bigBedToBed", str(path), str(out)], check=True)
    return out


def difference(ours, ref, mode, tmp_path):
    """Describe how our file differs from the reference.

    Parameters
    ----------
    ours, ref : pathlib.Path
        Our file and the reference file.
    mode : str
        `bytes`, `ref-cols-1-8`, `except-col-6`, `sorted-lists` or `bigbed`.
    tmp_path : pathlib.Path
        Scratch directory.

    Returns
    -------
    str or None
        First difference, None if the files match under `mode`.
    """
    if mode == "bigbed":
        ours = bigbed_text(ours, tmp_path / "ours.bed")
        ref = bigbed_text(ref, tmp_path / "ref.bed")
        mode = "bytes"
    if mode == "bytes":
        if same_bytes(ours, ref):
            return None
        return first_difference(ours, ref, strip, strip) or "same lines, other bytes"
    if mode == "ref-cols-1-8":
        return first_difference(ours, ref, strip, first_8_columns)
    if mode == "except-col-6":
        return first_difference(ours, ref, except_column_6, except_column_6)
    return first_difference(ours, ref, sort_lists, sort_lists)


def test_matches_reference(case, tmp_path):
    """Our output equals its reference file under the case's mode.

    Parameters
    ----------
    case : tuple of (pathlib.Path, pathlib.Path, str)
        Our file, the reference file and the comparison mode.
    tmp_path : pathlib.Path
        Scratch directory.
    """
    ours, ref, mode = case
    assert Path(ref).exists(), f"no reference file {ref}"
    detail = difference(ours, ref, mode, tmp_path)
    assert detail is None, detail
