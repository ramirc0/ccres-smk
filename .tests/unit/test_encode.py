"""Tests of the ENCODE file listing."""

import _encode


def test_inputs_list_every_file_once():
    """List each file once, skip `NA` and keep the extension of its use."""
    hotspots = [["ENCSR000AAA", "ENCFF000AAA", "ENCFF000AAB", "ENCFF000AAC", "1"]]
    marks = [
        ["ENCSR000AAA", "ENCFF000AAB", "cell_a"],
        ["ENCSR000AAD", "NA", "cell_b"],
        ["ENCSR000AAE", "ENCFF000AAE", "cell_b"],
    ]
    peaks = [
        ["ENCSR000AAF", "ENCFF000AAF", "Homo sapiens A549"],
        ["ENCSR000AAF", "ENCFF000AAF", "Homo sapiens A549", "ENCDO000AAA"],
    ]
    assert _encode.encode_inputs(hotspots, marks, peaks) == [
        ("ENCSR000AAA", "ENCFF000AAA", "bed.gz"),
        ("ENCSR000AAA", "ENCFF000AAB", "bigWig"),
        ("ENCSR000AAE", "ENCFF000AAE", "bigWig"),
        ("ENCSR000AAF", "ENCFF000AAF", "bed.gz"),
    ]
