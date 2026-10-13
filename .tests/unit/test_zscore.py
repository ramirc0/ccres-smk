"""Tests of the host-independent log and summation of _zscore."""

import os
import subprocess
import sys

import _zscore
import numpy as np
import pytest

REFERENCE = """
import math, sys, numpy as np
x = np.fromfile(sys.argv[1])
np.array([math.log(v, 10) for v in x.tolist()]).tofile(sys.argv[2])
"""


@pytest.mark.parametrize("seed", range(3))
def test_pairwise_sum_matches_numpy(seed):
    """pairwise_sum equals np.add.reduce bit for bit up to one 8192 buffer."""
    rng = np.random.default_rng(seed)
    sizes = [*range(1, 300), *rng.integers(300, 8193, 200).tolist(), 8192]
    for n in sizes:
        a = rng.normal(0, 1, n) * 10.0 ** rng.integers(-3, 4, n)
        assert _zscore.pairwise_sum(a) == np.add.reduce(a), n


def test_buffered_sum_chunks():
    """buffered_sum adds pairwise sums of consecutive 8192-element buffers."""
    a = np.random.default_rng(0).normal(0, 1, 3 * 8192 + 17)
    expected = -0.0
    for start in range(0, len(a), 8192):
        expected += float(np.add.reduce(a[start : start + 8192]))
    assert _zscore.buffered_sum(a) == expected


def test_log10_matches_glibc_without_fma(tmp_path):
    """log10 equals math.log(x, 10) under glibc's non-FMA log."""
    rng = np.random.default_rng(1)
    mantissa = rng.integers(100000, 1000000, 200_000)
    exponent = rng.integers(-11, 1, 200_000)
    x = np.array([float(f"{m}e{e}") for m, e in zip(mantissa, exponent)])
    x = np.concatenate([x, [1.0, 10.0, 0.9375, 1.0644]])
    x.tofile(tmp_path / "x.bin")
    subprocess.run(
        [sys.executable, "-c", REFERENCE, tmp_path / "x.bin", tmp_path / "y.bin"],
        env={**os.environ, "GLIBC_TUNABLES": "glibc.cpu.hwcaps=-FMA"},
        check=True,
    )
    expected = np.fromfile(tmp_path / "y.bin")
    assert np.array_equal(_zscore.log10(x).view(np.int64), expected.view(np.int64))
