"""Collect (ours, reference, mode) cases for the reference comparison.

Each step maps our files, found by glob under `<results>/<genome>/`, to files
under `<reference>/V7-<genome>/`. Reference files with no counterpart of ours
are ignored, since the reference tree also holds later analyses.
"""

from pathlib import Path

import pytest

REFERENCE = "/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry"

# step -> (glob under ours, reference path template, mode[, genomes]).
# Templates take {g} (genome), {name} (our basename) and {stem} (our basename
# without its last suffix).
STEPS = {
    # The release's lists were edited later. hg38's Filtered-Mod keeps the May
    # 2022 rows in columns 1-8. mm10's Filtered.txt differs only in column 6,
    # where a 2023 edit turned "," into ":".
    "tf_list": [
        (
            "TF/{g}-TF-List.Filtered.txt",
            "{g}-TF/{g}-TF-List.Filtered-Mod.txt",
            "ref-cols-1-8",
            {"hg38"},
        ),
        (
            "TF/{g}-TF-List.Filtered.txt",
            "{g}-TF/{g}-TF-List.Filtered.txt",
            "except-col-6",
            {"mm10"},
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
    # The release keeps no hg38 tmp.bed.
    "tf": [
        ("TF/tmp.bed", "{g}-TF/{name}", "bytes", {"mm10"}),
        ("TF/tmp.summits", "{g}-TF/{name}", "bytes"),
        ("TF/tmp.rPeaks", "{g}-TF/{name}", "bytes"),
        ("TF/{g}-TF-Clusters.txt", "{g}-TF/{name}", "sorted-lists"),
        ("TF/tmp.no-overlap", "{g}-TF/{name}", "sorted-lists"),
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
        ),
        (
            "Cell-Type-Specific/Individual-Files/*.bed.gz",
            "Cell-Type-Specific/Individual-Files/{stem}",
            "bytes",
        ),
        (
            "Cell-Type-Specific/Individual-Files/*.bigBed",
            "Cell-Type-Specific/Individual-Files/{name}",
            "bytes",
        ),
    ],
    # The release has no mm10 PLS gene list.
    "pls": [("{g}-PLS-Gene-List.txt", "PLS-Gene-List.{g}-V4.txt", "bytes", {"hg38"})],
}


# D5: the release names the neural crest cell's files without its H3K4me3 "NA"
# and has none for GM12866, whose only mark is that "NA".
RENAMED = {
    "hg38/Cell-Type-Specific/Individual-Files/ENCFF592CNZ_NA_ENCFF346DEM_ENCFF192FOK.bed.gz": (
        "Cell-Type-Specific/Individual-Files/ENCFF592CNZ_ENCFF346DEM_ENCFF192FOK.bed"
    ),
    "hg38/Cell-Type-Specific/Individual-Files/ENCFF592CNZ_NA_ENCFF346DEM_ENCFF192FOK.bigBed": (
        "Cell-Type-Specific/Individual-Files/ENCFF592CNZ_ENCFF346DEM_ENCFF192FOK.bigBed"
    ),
}
UNMATCHED = {
    "hg38/Cell-Type-Specific/Individual-Files/NA.bed.gz",
    "hg38/Cell-Type-Specific/Individual-Files/NA.bigBed",
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
            for pattern, template, mode, *genomes in rules:
                if genomes and g not in genomes[0]:
                    continue
                for path in sorted(ours_dir.glob(pattern.format(g=g))):
                    case_id = f"{g}/{path.relative_to(ours_dir)}"
                    ref = RENAMED.get(
                        case_id, template.format(g=g, name=path.name, stem=path.stem)
                    )
                    marks = [getattr(pytest.mark, step)]
                    if case_id in UNMATCHED:
                        marks.append(pytest.mark.skip(reason="D5: no release file"))
                    params.append(
                        pytest.param(
                            (path, reference / f"V7-{g}" / ref, mode),
                            marks=marks,
                            id=case_id,
                        )
                    )
    metafunc.parametrize("case", params)
