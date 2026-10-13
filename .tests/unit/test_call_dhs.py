"""Differential test of call_dhs against the loop of 0_Call-DHSs.sh."""

import subprocess

import call_dhs
import numpy as np
import pytest

MIN_P = 10
# Fixed-notation FDRs as Hotspot2 prints them. 0.0099999 rounds to order 2.
FIXED = ["0.0468394", "0.0099999", "0.00999", "0.000277824", "0"]


def track(rng, regions):
    """Draw a Hotspot2-like track: runs of short positions with FDRs.

    Parameters
    ----------
    rng : numpy.random.Generator
        Random source.
    regions : int
        Runs per chromosome.

    Returns
    -------
    str
        Track text, sorted by position within each chromosome.
    """
    lines = []
    for chrom in ["chr1", "chr2", "chr10"]:
        end = 0
        for _ in range(regions):
            end += int(rng.choice([1, 2, 30, 500]))
            # Runs whose floor exceeds MIN_P can stay over 350 bp to the end.
            floor = rng.integers(1, MIN_P + 3)
            for _ in range(rng.integers(1, 60)):
                start = end + int(rng.choice([0] * 8 + [1, 2]))
                end = start + int(rng.integers(1, 30))
                if rng.random() < 0.1:
                    fdr = str(rng.choice(FIXED))
                else:
                    order = floor + rng.integers(0, 4)
                    fdr = f"{rng.uniform(1, 10) * 10.0**-order:g}"
                lines.append(f"{chrom}\t{start}\t{end}\ti\t{fdr}\n")
    return "".join(lines)


def original_sweep(hotspots, tmp_path, bedtools, toolkit):
    """Run the threshold sweep of 0_Call-DHSs.sh up to 1E-MIN_P.

    Parameters
    ----------
    hotspots : pathlib.Path
        Track to call.
    tmp_path : pathlib.Path
        Working directory.
    bedtools : callable
        bedtools 2.30.0 runner.
    toolkit : pathlib.Path
        Original Toolkit directory.

    Returns
    -------
    tuple of str
        Text of the DHS and excluded files.
    """
    kept, merged = hotspots, {}
    for j in range(2, MIN_P + 1):
        filtered = subprocess.run(
            ["python3", toolkit / "filter-long-double.py", kept, f"1E-{j}"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        kept = tmp_path / f"tmp.{j}"
        kept.write_text(filtered)
        merged[j] = tmp_path / f"1E-{j}.bed"
        merged[j].write_text(
            "".join(
                line
                for line in bedtools(
                    "merge", "-d", 1, "-c", 5, "-o", "min", "-i", kept
                ).splitlines(keepends=True)
                if length(line) >= 50
            )
        )
    peaks = tmp_path / "peaks"
    peaks.write_text(short(merged[2].read_text()))
    for j in range(3, MIN_P + 1):
        new = bedtools("intersect", "-v", "-a", merged[j], "-b", peaks)
        with open(peaks, "a") as f:
            f.write(short(new))
    excluded = bedtools("intersect", "-v", "-a", merged[MIN_P], "-b", peaks)
    return peaks.read_text(), excluded


def length(line):
    """Return the length of a BED line.

    Parameters
    ----------
    line : str
        BED line.

    Returns
    -------
    int
        End minus start.
    """
    fields = line.split("\t")
    return int(fields[2]) - int(fields[1])


def short(text):
    """Keep the lines of at most 350 bp, as the original's awk filter.

    Parameters
    ----------
    text : str
        BED text.

    Returns
    -------
    str
        Lines of at most 350 bp.
    """
    return "".join(
        line for line in text.splitlines(keepends=True) if length(line) <= 350
    )


@pytest.mark.parametrize("seed", range(4))
def test_sweep(bedtools, original, tmp_path, seed):
    """call_dhs reproduces the DHS and excluded files of the original sweep."""
    hotspots = tmp_path / "hotspots.bed"
    hotspots.write_text(track(np.random.default_rng(seed), 100))
    dhs, excluded = tmp_path / "ours.DHSs.bed", tmp_path / "ours.Excluded.bed"
    call_dhs.main(
        [
            "--hotspots",
            str(hotspots),
            "--min-p",
            str(MIN_P),
            "--dhs",
            str(dhs),
            "--excluded",
            str(excluded),
        ]
    )
    expected_dhs, expected_excluded = original_sweep(
        hotspots, tmp_path, bedtools, original
    )
    assert expected_dhs and expected_excluded
    assert dhs.read_text() == expected_dhs
    assert excluded.read_text() == expected_excluded
