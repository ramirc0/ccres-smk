"""Check the frozen lists against the release's lists and the metadata snapshot.

Select with `-m inputs`. A frozen list keeps the release's rows in order, less
those whose file ENCODE withdrew, and names only released files.
"""

from functools import cache
from pathlib import Path

import polars as pl
import pytest

pytestmark = pytest.mark.inputs

GENOMES = ("hg38", "mm10")
# Frozen list -> the release's list under V7-<genome>/.
LISTS = {
    **{f"{m}-List.txt": f"{m}-List.txt" for m in ("DNase", "H3K4me3", "H3K27ac")},
    **{f"{m}-List.txt": f"{m}-List.txt" for m in ("CTCF", "ATAC")},
    "TF-List-All.txt": "{g}-TF/{g}-TF-List.All.txt",
}
HOTSPOTS = {
    "hg38": "hg38-rDHS/hg38-Hotspot-List.txt",
    "mm10": "mm10-rDHS/rDHS-cDHS-List.05-04-2022.mm10.txt",
}
FRIP = 4


def rows(path):
    """Read a tab-separated list.

    Parameters
    ----------
    path : pathlib.Path
        List file.

    Returns
    -------
    list of list of str
        Fields of each line.
    """
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f]


@cache
def files(meta):
    """Read the experiment and status of every file in the snapshot.

    Parameters
    ----------
    meta : str
        Directory of the ENCODE metadata snapshot.

    Returns
    -------
    dict of str to tuple of (str, str)
        Experiment and status per file accession.
    """
    table = pl.read_parquet(
        Path(meta) / "file.parquet", columns=["accession", "dataset", "status"]
    )
    exp = pl.col("dataset").str.extract(r"^/experiments/(\w+)/$")
    return {
        acc: (e, s) for acc, e, s in table.select("accession", exp, "status").rows()
    }


@pytest.fixture
def dirs(request):
    """Return the frozen input, release and metadata directories.

    Parameters
    ----------
    request : pytest.FixtureRequest
        Gives the command-line options.

    Returns
    -------
    tuple of (pathlib.Path, pathlib.Path, str)
        Frozen inputs, release trees and metadata snapshot.
    """
    option = request.config.getoption
    return Path(option("inputs")), Path(option("reference")), option("metadata")


def kept(row, name):
    """Return the columns of a list row that the refresh keeps.

    Parameters
    ----------
    row : list of str
        List row.
    name : str
        Frozen list name.

    Returns
    -------
    list of str
        Columns after the file, less a TF row's FRiP.
    """
    rest = row[2:]
    return rest[: FRIP - 2] + rest[FRIP - 1 :] if name == "TF-List-All.txt" else rest


@pytest.mark.parametrize("name", [*LISTS, "Hotspot-List.txt"])
@pytest.mark.parametrize("genome", GENOMES)
def test_lists_name_released_files(genome, name, dirs):
    """Name only released files, each under its experiment now."""
    inputs, _, meta = dirs
    found = files(meta)
    columns = slice(1, 4) if name == "Hotspot-List.txt" else slice(1, 2)
    bad = [
        (r[0], acc, found.get(acc))
        for r in rows(inputs / genome / name)
        for acc in r[columns]
        if acc != "NA" and found.get(acc) != (r[0], "released")
    ]
    assert not bad


@pytest.mark.parametrize("name", LISTS)
@pytest.mark.parametrize("genome", GENOMES)
def test_lists_keep_the_release_rows(genome, name, dirs):
    """Keep each release row in order, less those of withdrawn files."""
    inputs, reference, meta = dirs
    found = files(meta)
    frozen = rows(inputs / genome / name)
    i = 0
    for r in rows(reference / f"V7-{genome}" / LISTS[name].format(g=genome)):
        status = found.get(r[1], (None, None))[1]
        if r[1] != "NA" and status not in ("released", "archived"):
            continue
        assert kept(frozen[i], name) == kept(r, name), r
        assert frozen[i][1] == r[1] or status == "archived", r
        i += 1
    assert i == len(frozen)


@pytest.mark.parametrize("genome", GENOMES)
def test_hotspot_list_extends_the_release(genome, dirs):
    """Start with the release's rows and cover every DNase list experiment."""
    inputs, reference, meta = dirs
    found = files(meta)
    frozen = rows(inputs / genome / "Hotspot-List.txt")
    release = [
        r
        for r in rows(reference / f"V7-{genome}" / HOTSPOTS[genome])
        if found[r[1]][1] in ("released", "archived")
    ]
    for old, new in zip(release, frozen):
        assert new[1] == old[1] or found[old[1]][1] == "archived", old
    dnase = {r[0] for r in rows(inputs / genome / "DNase-List.txt")}
    assert dnase <= {r[0] for r in frozen}
