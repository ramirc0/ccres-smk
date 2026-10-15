"""Shared helpers for the unit tests: script imports and the external oracles."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[2] / "workflow" / "scripts"))


def run(args, text=None, env=None):
    """Run a command and return its stdout.

    Parameters
    ----------
    args : list of str
        Command and arguments.
    text : str, optional
        Standard input.
    env : dict of str to str, optional
        Variables added to the environment.

    Returns
    -------
    str
        Standard output.
    """
    return subprocess.run(
        args,
        input=text,
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, **(env or {})},
    ).stdout


@pytest.fixture
def bedtools():
    """Run bedtools 2.30.0, the version of the original run.

    Returns
    -------
    callable
        Takes bedtools arguments, returns its stdout.
    """
    binary = os.environ["BEDTOOLS"]
    version = run([binary, "--version"]).strip()
    assert version == "bedtools v2.30.0", version
    return lambda *args: run([binary, *map(str, args)])


@pytest.fixture
def gnu_sort():
    """Run GNU sort under the C locale, which compares bytes.

    Returns
    -------
    callable
        Takes sort arguments and the input text, returns sorted text.
    """
    return lambda args, text: run(["sort", *args], text, {"LC_ALL": "C"})


ORIGINAL = Path(
    os.environ.get(
        "CCRE_ORIGINAL", "~/Git_Repos/ENCODE-cCREs/Version-4/cCRE-Pipeline/Toolkit"
    )
).expanduser()


@pytest.fixture
def original():
    """Locate the original Toolkit scripts, skipping when they are absent.

    Returns
    -------
    pathlib.Path
        Toolkit directory of the original pipeline.
    """
    if not ORIGINAL.is_dir():
        pytest.skip(f"original Toolkit not found at {ORIGINAL}")
    return ORIGINAL
