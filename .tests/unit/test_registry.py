"""Differential tests of _registry against the original Toolkit scripts."""

import subprocess

import _bed
import _registry
import numpy as np
import polars as pl
import pytest

SIGNALS = ["1.5", "1.50", "0.75", "2", "0.3", "1.18826"]


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
        Columns chrom, start, end, name, fdr, signal in GNU sort order.
    """
    start = rng.integers(0, 30000, n)
    ids = rng.integers(0, n // 4, n) if duplicate_names else np.arange(n)
    df = pl.DataFrame(
        {
            "chrom": rng.choice(["chr1", "chr2", "chr10"], n).tolist(),
            "start": start.astype(str).tolist(),
            "end": (start + rng.integers(150, 351, n)).astype(str).tolist(),
            "name": [f"ENCFF{i % 7:03d}AAA-{i}" for i in ids],
            "fdr": ["1e-05"] * n,
            "signal": rng.choice(SIGNALS, n).tolist(),
        }
    )
    return df[_bed.gnu_sort(df, text_keys=["chrom"], numeric_keys=["start"])]


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
    _bed.write_tsv(df, sorted_path)
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
    codes = _bed.chrom_codes(df["chrom"])
    picks = _registry.representative_peaks(
        codes,
        df["start"].cast(pl.Int64).to_numpy(),
        df["end"].cast(pl.Int64).to_numpy(),
        np.array([float(s) for s in df["signal"]]),
        df["name"],
    )
    ours = "".join(f"{line}\n" for line in _bed.lines(df[picks]).to_list())
    assert ours == original_loop(df, tmp_path, bedtools, original)


def test_accessions(original, tmp_path):
    """Accessions reproduces make-region-accession.py without key collisions."""
    rng = np.random.default_rng(0)
    start = rng.integers(1, 5000, 400)
    regions = [("chr1", str(s), str(s + 200)) for s in start]
    regions += [("chr1", "10033", "10250")]
    previous = [
        ("chr1", str(s), str(s + 200), f"EH37D{i:07d}", "V1")
        for i, s in enumerate(start[::3])
    ]
    previous += [
        ("chr1", "10033", "10250", "EH37D0999999", "V2"),
        ("chr1", "12", "3", "EH37D0000007", "V3"),
    ]
    previous += [previous[0][:3] + ("EH37D0000123", "V3")]
    region_path, previous_path = tmp_path / "regions.bed", tmp_path / "previous.bed"
    region_path.write_text("".join("\t".join(r) + "\n" for r in regions))
    previous_path.write_text("".join("\t".join(r) + "\n" for r in previous))

    known = _registry.previous_accessions(previous, "rDHS")
    ours = _registry.accessions(regions, "hg38", "rDHS", known)
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
    assert "".join("\t".join(r) + f"\t{a}\n" for r, a in zip(regions, ours)) == expected


def test_accessions_key_on_coordinates():
    """Reuse a previous accession only for the same chrom, start and end."""
    previous = [
        ("chr1", "10033", "10250", "EH37D0999999", "V2"),
        ("chr1", "12", "3", "EH37D0000007", "V3"),
    ]
    regions = [("chr1", "1003", "310250"), ("chr11", "2", "3"), ("chr1", "12", "3")]
    known = _registry.previous_accessions(previous, "rDHS")
    assert _registry.accessions(regions, "hg38", "rDHS", known) == [
        "EH38D1000000",
        "EH38D1000001",
        "EH38D0000007",
    ]
