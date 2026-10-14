"""Locate ENCODE files: the local mirror first, then the download directory."""

from pathlib import Path


def encode_path(encode_root, download_dir, exp, acc, ext):
    """Return the path of an ENCODE file.

    Parameters
    ----------
    encode_root : str
        Local mirror, laid out as `<exp>/<acc>.<ext>`.
    download_dir : str
        Same layout, for files downloaded because the mirror lacks them.
    exp : str
        Experiment accession (`ENCSR...`).
    acc : str
        File accession (`ENCFF...`).
    ext : str
        File extension, such as `bed.gz` or `bigWig`.

    Returns
    -------
    str
        The mirror path if it exists, else the download path.
    """
    local = Path(encode_root, exp, f"{acc}.{ext}")
    return str(local) if local.exists() else f"{download_dir}/{exp}/{acc}.{ext}"


def encode_inputs(hotspots, marks, peaks):
    """List the ENCODE files that list rows name, each once.

    Parameters
    ----------
    hotspots : iterable of list of str
        Hotspot list rows: experiment, peak file, signal file, then others.
    marks : iterable of list of str
        Mark list rows: experiment, then signal file. `NA` names no file.
    peaks : iterable of list of str
        TF list rows: experiment, then peak file.

    Returns
    -------
    list of tuple of str
        Experiment, file accession and extension of each file, sorted.
    """
    files = set()
    for exp, peak, signal, *_ in hotspots:
        files.update([(exp, peak, "bed.gz"), (exp, signal, "bigWig")])
    files.update((r[0], r[1], "bigWig") for r in marks if r[1] != "NA")
    files.update((r[0], r[1], "bed.gz") for r in peaks)
    return sorted(files)
