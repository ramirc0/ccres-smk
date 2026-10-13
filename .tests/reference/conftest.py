"""Collect (ours, reference, mode) cases for the reference comparison.

Each step maps our files, found by glob under `<results>/<genome>/`, to files
under `<reference>/V7-<genome>/`. Reference files with no counterpart of ours
are ignored, since the reference tree also holds later analyses.
"""

from pathlib import Path

import pytest

REFERENCE = "/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry"

# step -> (glob under ours, reference path template, mode).
# Templates take {g} (genome), {name} (our basename) and {stem} (our basename
# without its last suffix).
STEPS = {
    "tf_list": [
        (
            "TF/{g}-TF-List.Filtered.txt",
            "{g}-TF/{g}-TF-List.Filtered-Mod.txt",
            "ref-cols-1-8",
        ),
    ],
    "dhs": [
        ("rDHS/DHSs/*.bed", "{g}-rDHS/DHSs/{name}", "bytes"),
        ("rDHS/Processed-DHSs/*.bed", "{g}-rDHS/Processed-DHSs/output.{stem}", "bytes"),
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
    "tf": [
        ("TF/tmp.bed", "{g}-TF/{name}", "bytes"),
        ("TF/tmp.summits", "{g}-TF/{name}", "bytes"),
        ("TF/tmp.rPeaks", "{g}-TF/{name}", "bytes"),
        ("TF/{g}-TF-Clusters.txt", "{g}-TF/{name}", "sorted-lists"),
        ("TF/tmp.no-overlap", "{g}-TF/{name}", "sorted-lists"),
        ("TF/{g}-rTFBS-Summary.txt", "{g}-TF/{name}", "sorted-lists"),
        ("TF/{g}-rTFBS.bed", "{g}-TF/{name}", "bytes"),
    ],
    "anchors": [("{g}-Anchors.bed", "{name}", "bytes")],
    "signal": [("signal-output/*/*.txt", "signal-output/{name}", "bytes")],
    "maxz": [("{g}-*-maxZ.txt", "{name}", "bytes")],
    "ccres": [("{g}-cCREs-Unfiltered.bed", "{name}", "bytes")],
    "cts": [
        (
            "Cell-Type-Specific/Master-Cell-List.txt",
            "Cell-Type-Specific/{name}",
            "bytes",
        ),
        (
            "Cell-Type-Specific/Individual-Files/*.bed",
            "Cell-Type-Specific/Individual-Files/{name}",
            "bytes",
        ),
        (
            "Cell-Type-Specific/Individual-Files/*.bigBed",
            "Cell-Type-Specific/Individual-Files/{name}",
            "bigbed",
        ),
    ],
    "pls": [("{g}-PLS-Gene-List.txt", "PLS-Gene-List.{g}-V4.txt", "bytes")],
}


def pytest_addoption(parser):
    """Add the result and reference directory options.

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


def pytest_configure(config):
    """Register one marker per step, for selection with `-m <step>`.

    Parameters
    ----------
    config : pytest.Config
        Test session config.
    """
    for step in STEPS:
        config.addinivalue_line("markers", f"{step}: outputs of step {step}")


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
            for pattern, template, mode in rules:
                for path in sorted(ours_dir.glob(pattern.format(g=g))):
                    ref = template.format(g=g, name=path.name, stem=path.stem)
                    params.append(
                        pytest.param(
                            (path, reference / f"V7-{g}" / ref, mode),
                            marks=getattr(pytest.mark, step),
                            id=f"{g}/{path.relative_to(ours_dir)}",
                        )
                    )
    metafunc.parametrize("case", params)
