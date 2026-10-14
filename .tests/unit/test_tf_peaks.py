"""Tests of the TF peak resizing of tf_peaks."""

import polars as pl
import tf_peaks


def test_resized_peaks_keep_their_summit():
    """A resized peak's offset points at the summit of the original peak."""
    peaks = pl.DataFrame(
        {
            "c1": ["chr1"] * 3,
            "c2": ["1000", "5000", "9000"],
            "c3": ["1100", "5500", "9200"],
            "c4": ["."] * 3,
            "c5": ["0"] * 3,
            "c6": ["."] * 3,
            "c7": ["12.5"] * 3,
            "c8": ["-1"] * 3,
            "c9": ["4.2"] * 3,
            "c10": ["30", "400", "50"],
        }
    )
    bed = tf_peaks.resized(peaks, "ENCSR000AAA", "CTCF", "K562")
    summits = tf_peaks.summit_points(bed)
    assert bed["c10"].to_list() == ["75", "175", "50"]
    assert summits["c2"].to_list() == ["1030", "5400", "9050"]
    assert summits["c3"].to_list() == ["1031", "5401", "9051"]
