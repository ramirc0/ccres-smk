"""Collect (ours, reference, mode) cases for the reference comparison.

Each step maps our files, found by glob under `<results>/<genome>/`, to files
under `<reference>/V7-<genome>/`. Reference files with no counterpart of ours
are ignored, since the reference tree also holds later analyses.
"""

import re
from functools import cache
from pathlib import Path

import pytest

REFERENCE = "/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry"
METADATA = "/home/ramirezc/Projects/encode-metadata-smk/results/2026-10-09/parquets"

# step -> (glob under ours, reference path template, mode[, genomes]).
# Templates take {g} (genome), {name} (the release's basename of our file) and
# {stem} (that name without its last suffix).
STEPS = {
    # The release's lists were edited later. hg38's Filtered-Mod keeps the May
    # 2022 rows in columns 1-8. mm10's Filtered.txt differs only in column 6,
    # where a 2023 edit turned "," into ":".
    "tf_list": [
        (
            "TF/{g}-TF-List-Filtered.txt",
            "{g}-TF/{g}-TF-List.Filtered-Mod.txt",
            "ref-cols-1-8",
            {"hg38"},
        ),
        (
            "TF/{g}-TF-List-Filtered.txt",
            "{g}-TF/{g}-TF-List.Filtered.txt",
            "except-col-6",
            {"mm10"},
        ),
    ],
    "dhs": [
        ("rDHS/DHSs/*.bed", "{g}-rDHS/DHSs/{name}", "bytes"),
        ("rDHS/Processed-DHSs/*.bed", "{g}-rDHS/Processed-DHSs/{name}", "bytes"),
    ],
    "rdhs": [
        ("rDHS/{g}-DHS-All.bed", "{g}-rDHS/{name}", "bytes"),
        ("rDHS/{g}-DHS-Filtered.bed", "{g}-rDHS/{name}", "bytes"),
        ("rDHS/{g}-rDHS-PreStam.txt", "{g}-rDHS/{name}", "bytes"),
        ("rDHS/completeCoverage.bed", "{g}-rDHS/{name}", "bytes"),
        ("rDHS/NoCoverage.bed", "{g}-rDHS/{name}", "bytes"),
        ("rDHS/{g}-rDHS-Filtered-Summary.txt", "{g}-rDHS/{name}", "bytes"),
        ("rDHS/{g}-rDHS-Filtered.bed", "{g}-rDHS/{name}", "bytes"),
    ],
    # The release keeps no hg38 tmp.bed.
    "tf": [
        ("TF/tmp.bed", "{g}-TF/{name}", "bytes", {"mm10"}),
        ("TF/tmp-summits", "{g}-TF/{name}", "bytes"),
        ("TF/tmp-rPeaks", "{g}-TF/{name}", "bytes"),
        ("TF/{g}-TF-Clusters.txt", "{g}-TF/{name}", "sorted-lists"),
        ("TF/tmp-no-overlap", "{g}-TF/{name}", "sorted-lists"),
        ("TF/{g}-rTFBS-Summary.txt", "{g}-TF/{name}", "sorted-lists"),
        ("TF/{g}-rTFBS.bed", "{g}-TF/{name}", "bytes"),
    ],
    "anchors": [("{g}-Anchors.bed", "{name}", "bytes")],
    # The release's ATAC files come from a September 2023 run whose glibc used
    # its FMA log. The port uses the non-FMA log of the May 2022 run (D4).
    "signal": [
        *(
            (f"signal-output/{mark}/*.txt", "signal-output/{name}", "bytes")
            for mark in ("DNase", "H3K4me3", "H3K27ac", "CTCF")
        ),
        ("signal-output/ATAC/*.txt", "signal-output/{name}", "z-close"),
    ],
    "maxz": [
        *(
            (f"{{g}}-{mark}-maxZ.txt", "{name}", "bytes")
            for mark in ("DNase", "H3K4me3", "H3K27ac", "CTCF")
        ),
        ("{g}-ATAC-maxZ.txt", "{name}", "z-close"),
    ],
    "ccres": [("{g}-cCREs-Unfiltered.bed", "{name}", "bytes")],
    "cts": [
        (
            "Cell-Type-Specific/Master-Cell-List.txt",
            "Cell-Type-Specific/{name}",
            "bytes",
            {"hg38"},
        ),
        (
            "Cell-Type-Specific/Individual-Files/*.bed.gz",
            "Cell-Type-Specific/Individual-Files/{stem}",
            "bytes",
            {"hg38"},
        ),
        (
            "Cell-Type-Specific/Individual-Files/*.bigBed",
            "Cell-Type-Specific/Individual-Files/{name}",
            "bytes",
            {"hg38"},
        ),
        # mm10's step 7 was rebuilt in 2025. Old/ keeps the May 2022 files, with
        # each BED gzipped in November 2022.
        (
            "Cell-Type-Specific/Master-Cell-List.txt",
            "Cell-Type-Specific/Old/{name}",
            "bytes",
            {"mm10"},
        ),
        (
            "Cell-Type-Specific/Individual-Files/*.bed.gz",
            "Cell-Type-Specific/Old/Individual-Files/{name}",
            "bytes",
            {"mm10"},
        ),
        (
            "Cell-Type-Specific/Individual-Files/*.bigBed",
            "Cell-Type-Specific/Old/Individual-Files/{name}",
            "bytes",
            {"mm10"},
        ),
    ],
    # The release has no mm10 PLS gene list.
    "pls": [("{g}-PLS-Gene-List.txt", "PLS-Gene-List.{g}-V4.txt", "bytes", {"hg38"})],
}


# Our basenames that differ from the release's, apart from step 7's.
RELEASE_NAMES = [
    (re.compile(r"^(ENCFF\w{6})-(DHSs|Excluded)\.bed$"), r"\1.\2.bed"),
    (re.compile(r"^(ENCFF\w{6})-(\d+)\.bed$"), r"output.\1.\2"),
    (re.compile(r"^tmp-"), "tmp."),
]


@cache
def split_names(ours_dir):
    """Map each step 7 file stem of a run to the name the Split scripts gave it.

    Parameters
    ----------
    ours_dir : pathlib.Path
        Run directory of one genome.

    Returns
    -------
    dict of str to str
        Row ID, its first experiment, to the row's file accessions joined by
        `_`.
    """
    with open(ours_dir / "Cell-Type-Specific" / "Master-Cell-List.txt") as f:
        rows = [line.split("\t") for line in f]
    return {
        next(e for e in r[0:8:2] if e != "---"): "_".join(
            a for a in r[1:8:2] if a != "---"
        )
        for r in rows
    }


def release_name(path, ours_dir):
    """Return the release's basename of one of our files.

    Parameters
    ----------
    path : pathlib.Path
        Our file.
    ours_dir : pathlib.Path
        Run directory of its genome.

    Returns
    -------
    str
        The basename the release gave the same output.
    """
    if path.parent.name == "Individual-Files":
        row_id, suffix = path.name.split(".", 1)
        return f"{split_names(ours_dir)[row_id]}.{suffix}"
    name = path.name
    for pattern, replacement in RELEASE_NAMES:
        name = pattern.sub(replacement, name)
    return name


def pytest_addoption(parser):
    """Add the result, reference, frozen input and metadata directory options.

    Parameters
    ----------
    parser : pytest.Parser
        Option parser.
    """
    parser.addoption(
        "--results",
        default="data/results/v7",
        help="Run directory holding one subdirectory per genome.",
    )
    parser.addoption(
        "--reference",
        default=REFERENCE,
        help="Directory holding the V7-<genome> reference trees.",
    )
    parser.addoption(
        "--inputs",
        default="data/inputs-2026-10-09",
        help="Frozen inputs, one subdirectory per genome.",
    )
    parser.addoption(
        "--metadata", default=METADATA, help="ENCODE metadata snapshot of the inputs."
    )


def pytest_configure(config):
    """Register one marker per step and one for the inputs, for `-m <marker>`.

    Parameters
    ----------
    config : pytest.Config
        Test session config.
    """
    for step in STEPS:
        config.addinivalue_line("markers", f"{step}: outputs of step {step}")
    config.addinivalue_line("markers", "inputs: frozen lists against the release's")


def pytest_generate_tests(metafunc):
    """Parametrize `case` over every output of every step and genome.

    Parameters
    ----------
    metafunc : pytest.Metafunc
        Collection hook argument.
    """
    if "case" not in metafunc.fixturenames:
        return
    results = Path(metafunc.config.getoption("results"))
    reference = Path(metafunc.config.getoption("reference"))
    params = []
    for ours_dir in sorted(p for p in results.glob("*") if p.is_dir()):
        g = ours_dir.name
        for step, rules in STEPS.items():
            for pattern, template, mode, *genomes in rules:
                if genomes and g not in genomes[0]:
                    continue
                for path in sorted(ours_dir.glob(pattern.format(g=g))):
                    name = release_name(path, ours_dir)
                    ref = template.format(g=g, name=name, stem=name.rsplit(".", 1)[0])
                    params.append(
                        pytest.param(
                            (path, reference / f"V7-{g}" / ref, mode),
                            marks=getattr(pytest.mark, step),
                            id=f"{g}/{path.relative_to(ours_dir)}",
                        )
                    )
    metafunc.parametrize("case", params)
