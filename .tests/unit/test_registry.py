"""Differential tests of _registry against the original Toolkit scripts."""

import subprocess

import _bed
import _registry
import numpy as np
import polars as pl
import pytest

SIGNALS = [1.5, 0.75, 2.0, 0.3, 1.18826]
REGION = {"chrom": pl.String, "start": pl.Int64, "end": pl.Int64}
PREVIOUS = {**REGION, "accession": pl.String, "version": pl.String}


def peaks(rng, n, duplicate_names=False):
    """Draw sorted DHS-like peaks with signal ties.

    Parameters
    ----------
    rng : numpy.random.Generator
        Random source.
    n : int
        Number of peaks.
    duplicate_names : bool
        Reuse names, as an experiment listed twice would.

    Returns
    -------
    polars.DataFrame
        Columns chrom, start, end, name, fdr, signal, sorted.
    """
    start = rng.integers(0, 30000, n)
    ids = rng.integers(0, n // 4, n) if duplicate_names else np.arange(n)
    df = pl.DataFrame(
        {
            "chrom": rng.choice(["chr1", "chr2", "chr10"], n).tolist(),
            "start": start,
            "end": start + rng.integers(150, 351, n),
            "name": [f"ENCFF{i % 7:03d}AAA-{i}" for i in ids],
            "fdr": [1e-05] * n,
            "signal": rng.choice(SIGNALS, n),
        }
    )
    return df.sort("chrom", "start", "end", "name")


def original_loop(df, tmp_path, bedtools, toolkit):
    """Run the merge, pick and intersect loop of 2_Create-rDHSs.sh.

    Parameters
    ----------
    df : polars.DataFrame
        Sorted peaks.
    tmp_path : pathlib.Path
        Working directory.
    bedtools : callable
        bedtools 2.30.0 runner.
    toolkit : pathlib.Path
        Original Toolkit directory.

    Returns
    -------
    str
        Text of the rPeaks file.
    """
    sorted_path, merge, rpeaks = (
        tmp_path / "sorted",
        tmp_path / "merge",
        tmp_path / "rPeaks",
    )
    _bed.write_bed(df, sorted_path)
    rpeaks.write_text("")
    while sorted_path.stat().st_size:
        merge.write_text(
            bedtools("merge", "-i", sorted_path, "-c", "4,6", "-o", "collapse,collapse")
        )
        picks = subprocess.run(
            ["python3", toolkit / "pick-best-peak.py", merge],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        (tmp_path / "peak-list").write_text(picks)
        chosen = subprocess.run(
            [
                "awk",
                "FNR==NR {x[$1];next} ($4 in x)",
                tmp_path / "peak-list",
                sorted_path,
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        with open(rpeaks, "a") as f:
            f.write(chosen)
        sorted_path.write_text(
            bedtools("intersect", "-v", "-a", sorted_path, "-b", rpeaks)
        )
    return rpeaks.read_text()


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("duplicate_names", [False, True])
def test_representative_peaks(bedtools, original, tmp_path, seed, duplicate_names):
    """representative_peaks reproduces the original loop, ties and shared names included."""
    df = peaks(np.random.default_rng(seed), 1500, duplicate_names)
    picks = _registry.representative_peaks(df, "signal", "name")
    ours = picks.write_csv(separator="\t", include_header=False)
    assert ours == original_loop(df, tmp_path, bedtools, original)


def test_accessions(original, tmp_path):
    """Accessions reproduces make-region-accession.py without key collisions."""
    rng = np.random.default_rng(0)
    start = rng.integers(1, 5000, 400)
    regions = [("chr1", int(s), int(s) + 200) for s in start]
    regions += [("chr1", 10033, 10250)]
    previous = [
        ("chr1", int(s), int(s) + 200, f"EH37D{i:07d}", "V1")
        for i, s in enumerate(start[::3])
    ]
    previous += [
        ("chr1", 10033, 10250, "EH37D0999999", "V2"),
        ("chr1", 12, 3, "EH37D0000007", "V3"),
    ]
    previous += [previous[0][:3] + ("EH37D0000123", "V3")]
    regions = pl.DataFrame(regions, schema=REGION, orient="row")
    previous = pl.DataFrame(previous, schema=PREVIOUS, orient="row")
    region_path, previous_path = tmp_path / "regions.bed", tmp_path / "previous.bed"
    _bed.write_bed(regions, region_path)
    _bed.write_bed(previous, previous_path)

    ours = _registry.accessions(regions, "hg38", "rDHS", previous)
    expected = subprocess.run(
        [
            "python3",
            original / "make-region-accession.py",
            region_path,
            "hg38",
            "rDHS",
            previous_path,
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    ours = regions.with_columns(ours).write_csv(separator="\t", include_header=False)
    assert ours == expected


def test_accessions_key_on_coordinates():
    """Reuse a previous accession only for the same chrom, start and end."""
    previous = pl.DataFrame(
        [
            ("chr1", 10033, 10250, "EH37D0999999", "V2"),
            ("chr1", 12, 3, "EH37D0000007", "V3"),
        ],
        schema=PREVIOUS,
        orient="row",
    )
    regions = pl.DataFrame(
        [("chr1", 1003, 310250), ("chr11", 2, 3), ("chr1", 12, 3)],
        schema=REGION,
        orient="row",
    )
    assert _registry.accessions(regions, "hg38", "rDHS", previous).to_list() == [
        "EH38D1000000",
        "EH38D1000001",
        "EH38D0000007",
    ]
