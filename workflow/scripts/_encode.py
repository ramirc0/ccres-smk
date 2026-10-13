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
