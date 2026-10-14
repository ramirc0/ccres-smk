"""Tests of the step 3 TF experiment filter."""

import pytest
import tf_list


@pytest.mark.parametrize(
    ("frip", "peak", "expected"),
    [
        ("0.004", "ENCFF000AAA", True),
        ("1e-2", "ENCFF000AAA", True),
        ("0.003", "ENCFF000AAA", False),
        ("8.269033379039276e-05", "ENCFF000AAA", False),
        ("NA", "ENCFF000AAA", False),
        ("0.5", "NA", False),
    ],
)
def test_frip_must_be_numeric(frip, peak, expected):
    """Keep an experiment only if its FRiP is a number above 0.003."""
    fields = ["ENCSR000AAA", peak, "Homo sapiens K562", "CTCF", frip, "lab", "x"]
    assert tf_list.passes(fields) is expected
