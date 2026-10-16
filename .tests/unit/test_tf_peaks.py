"""Tests of the TF peak resizing of tf_peaks."""

import polars as pl
import tf_peaks


def test_resized_peaks_keep_their_summit():
    """A resized peak centers on the summit of the original peak."""
    peaks = pl.DataFrame(
        {
            "chrom": ["chr1"] * 3,
            "start": [1000, 5000, 9000],
            "end": [1100, 5500, 9200],
            "signal": [12.5] * 3,
            "offset": [30, 400, 50],
        }
    )
    bed = tf_peaks.resized(peaks, "ENCSR000AAA", "CTCF", "K562")
    assert bed["summit"].to_list() == [1030, 5400, 9050]
    assert bed["start"].to_list() == [955, 5225, 9000]
    assert bed["end"].to_list() == [1105, 5575, 9200]
    assert bed["name"].to_list() == [f"ENCSR000AAA-{i}" for i in (1, 2, 3)]
