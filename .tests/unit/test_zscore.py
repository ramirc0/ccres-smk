"""Tests of the host-independent log and the exact mean of _zscore."""

import os
import subprocess
import sys

import _zscore
import numpy as np

REFERENCE = """
import math, sys, numpy as np
x = np.fromfile(sys.argv[1])
np.array([math.log(v, 10) for v in x.tolist()]).tofile(sys.argv[2])
"""


def test_mean_std_sums_exactly():
    """Round each sum once, whatever the order of the values."""
    assert _zscore.mean_std(np.array([1e16, 1.0, -1e16, 3.0]))[0] == 1.0
    rng = np.random.default_rng(0)
    a = rng.normal(0, 1, 100_003) * 10.0 ** rng.integers(-3, 4, 100_003)
    assert _zscore.mean_std(a) == _zscore.mean_std(rng.permutation(a))


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
