"""Tests of the list refresh and the hotspot list that freeze_input derives."""

import freeze_input
import polars as pl
import pytest

FILE = {
    "accession": None,
    "dataset": None,
    "status": "released",
    "output_type": None,
    "file_format": None,
    "biological_replicates": "1",
    "genome_assembly": "GRCh38",
    "analyses": "/analyses/ENCAN000NEW/",
}


def metadata(tmp_path, files, experiments):
    """Write the file and experiment tables of a metadata snapshot.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to write the Parquet files in.
    files : list of dict
        File fields that differ from `FILE`.
    experiments : list of tuple of str
        Accession, default analysis, lab and biosample summary.

    Returns
    -------
    pathlib.Path
        The snapshot directory.
    """
    pl.DataFrame(
        [{**FILE, **f} for f in files], schema=dict.fromkeys(FILE, pl.String)
    ).write_parquet(tmp_path / "file.parquet")
    pl.DataFrame(
        experiments,
        schema=["accession", "default_analysis", "lab", "biosample_summary"],
        orient="row",
    ).write_parquet(tmp_path / "experiment.parquet")
    return tmp_path


def signal(accession, exp, **fields):
    """Return the fields of a signal bigWig.

    Parameters
    ----------
    accession, exp : str
        File and experiment accessions.
    **fields
        Fields that differ.

    Returns
    -------
    dict
        File fields.
    """
    return {
        "accession": accession,
        "dataset": f"/experiments/{exp}/",
        "output_type": "read-depth normalized signal",
        "file_format": "bigWig",
        **fields,
    }


def test_refresh_swaps_archived_and_drops_revoked(tmp_path):
    """Drop revoked files, swap archived ones and name each file's experiment."""
    old = {"status": "archived", "analyses": "/analyses/ENCAN000OLD/"}
    meta = metadata(
        tmp_path,
        [
            signal("ENCFF000AAA", "ENCSR000AAA", biological_replicates="1,2"),
            signal("ENCFF000AAB", "ENCSR000AAA", **old),
            signal("ENCFF000AAC", "ENCSR000AAA", biological_replicates="2"),
            signal("ENCFF000AAD", "ENCSR000AAA"),
            signal("ENCFF000AAE", "ENCSR000AAE", status="revoked"),
            signal("ENCFF000AAF", "ENCSR000AAF"),
        ],
        [("ENCSR000AAA", "/analyses/ENCAN000NEW/", "lab", "cell")]
        + [
            (e, "/analyses/ENCAN000NEW/", "lab", "cell")
            for e in ("ENCSR000AAE", "ENCSR000AAF")
        ],
    )
    rows = [
        ["ENCSR000AAA", "ENCFF000AAA", "cell_a"],
        ["ENCSR000AAA", "ENCFF000AAB", "cell_b"],
        ["ENCSR000AAE", "ENCFF000AAE", "cell_c"],
        ["ENCSR000OLD", "ENCFF000AAF", "cell_d"],
    ]
    assert freeze_input.refresh(rows, meta, "GRCh38") == [
        ["ENCSR000AAA", "ENCFF000AAA", "cell_a"],
        ["ENCSR000AAA", "ENCFF000AAD", "cell_b"],
        ["ENCSR000AAF", "ENCFF000AAF", "cell_d"],
    ]


def peaks(accession, output, replicates, **fields):
    """Return the fields of a peak file of ENCSR000AAA.

    Parameters
    ----------
    accession : str
        File accession.
    output : str
        Output type.
    replicates : str
        Biological replicates, comma-separated.
    **fields
        Fields that differ.

    Returns
    -------
    dict
        File fields.
    """
    return {
        "accession": accession,
        "dataset": "/experiments/ENCSR000AAA/",
        "output_type": output,
        "file_format": "bed",
        "biological_replicates": replicates,
        **fields,
    }


def test_refresh_takes_the_largest_replicate_subset(tmp_path):
    """Swap pooled optimal IDR peaks for the IDR peaks of one replicate."""
    old = {"status": "archived", "analyses": "/analyses/ENCAN000OLD/"}
    meta = metadata(
        tmp_path,
        [
            peaks("ENCFF000OPT", "optimal IDR thresholded peaks", "1,2", **old),
            peaks("ENCFF000ONE", "IDR thresholded peaks", "1"),
            peaks("ENCFF000CON", "conservative IDR thresholded peaks", "1"),
            peaks("ENCFF000OUT", "IDR thresholded peaks", "1,3"),
        ],
        [("ENCSR000AAA", "/analyses/ENCAN000NEW/", "lab", "cell")],
    )
    rows = [["ENCSR000AAA", "ENCFF000OPT", "cell_a"]]
    assert freeze_input.refresh(rows, meta, "GRCh38") == [
        ["ENCSR000AAA", "ENCFF000ONE", "cell_a"]
    ]


def test_refresh_rejects_two_replacements(tmp_path):
    """Fail when two files tie for an archived file's replicates."""
    old = {"status": "archived", "analyses": "/analyses/ENCAN000OLD/"}
    meta = metadata(
        tmp_path,
        [
            peaks("ENCFF000OPT", "IDR thresholded peaks", "1,2", **old),
            peaks("ENCFF000ONE", "IDR thresholded peaks", "1"),
            peaks("ENCFF000TWO", "IDR thresholded peaks", "2"),
        ],
        [("ENCSR000AAA", "/analyses/ENCAN000NEW/", "lab", "cell")],
    )
    with pytest.raises(ValueError, match="ENCFF000OPT"):
        freeze_input.refresh([["ENCSR000AAA", "ENCFF000OPT"]], meta, "GRCh38")


def test_tf_rows_take_the_frip_of_their_peak_file(tmp_path):
    """Write ENCODE's FRiP of a row's peak file, where ENCODE has one."""
    pl.DataFrame(
        {
            "files": ["/files/ENCFF000AAA/", "/files/ENCFF000AAB/"],
            "status": ["released", "deleted"],
            "frip_score_of_peaks": [0.015318143432016132, 0.5],
        }
    ).write_parquet(tmp_path / "chip_peak_enrichment_quality_metric.parquet")
    rows = [
        ["ENCSR000AAA", "ENCFF000AAA", "cell", "CTCF", "0.2", "lab"],
        ["ENCSR000AAB", "ENCFF000AAB", "cell", "CTCF", "0.1", "lab"],
        ["ENCSR000AAC", "NA", "cell", "CTCF", "NA", "lab"],
    ]
    assert freeze_input.encode_frip(rows, tmp_path) == [
        ["ENCSR000AAA", "ENCFF000AAA", "cell", "CTCF", "0.015318143432016132", "lab"],
        ["ENCSR000AAB", "ENCFF000AAB", "cell", "CTCF", "0.1", "lab"],
        ["ENCSR000AAC", "NA", "cell", "CTCF", "NA", "lab"],
    ]


def test_hotspot_rows_from_metadata(tmp_path):
    """Rebuild each listed replicate, then append the other replicates."""
    kinds = {"fdr": ("FDR cut rate", "bed"), "peaks": ("peaks", "starch")}

    def replicate(exp, rep, n):
        found = [signal(f"ENCFF{n}SIG", exp, biological_replicates=rep)]
        for kind, (output, form) in kinds.items():
            found.append(
                {
                    "accession": f"ENCFF{n}{kind[:3].upper()}",
                    "dataset": f"/experiments/{exp}/",
                    "output_type": output,
                    "file_format": form,
                    "biological_replicates": rep,
                }
            )
        return found

    bigbed = {"accession": "ENCFF9BBFDR", "dataset": "/experiments/ENCSR000AAA/"}
    meta = metadata(
        tmp_path,
        [
            *replicate("ENCSR000AAA", "1", 1),
            *replicate("ENCSR000AAA", "2", 2),
            *replicate("ENCSR000AAB", "1", 3),
            {**bigbed, "output_type": "FDR cut rate", "file_format": "bigBed"},
            {
                "accession": "ENCFF4FDR",
                "dataset": "/experiments/ENCSR000AAC/",
                "status": "revoked",
                "output_type": "FDR cut rate",
                "file_format": "bed",
            },
        ],
        [
            (e, "/analyses/ENCAN000NEW/", f"lab {e[-1]}", f"cell {e[-1]}")
            for e in ("ENCSR000AAA", "ENCSR000AAB", "ENCSR000AAC")
        ],
    )
    hotspots = [
        ["ENCSR000AAA", "ENCFF1FDR", "ENCFF1SIG", "ENCFF1PEA", "1", "x", "y"],
        ["ENCSR000AAC", "ENCFF4FDR", "ENCFF4SIG", "ENCFF4PEA", "1", "x", "y"],
    ]
    dnase = [["ENCSR000AAA", "ENCFF1SIG", "a"], ["ENCSR000AAB", "ENCFF3SIG", "b"]]
    assert freeze_input.hotspot_rows(hotspots, dnase, meta, "GRCh38") == [
        ["ENCSR000AAA", "ENCFF1FDR", "ENCFF1SIG", "ENCFF1PEA", "1", "lab A", "cell A"],
        ["ENCSR000AAA", "ENCFF2FDR", "ENCFF2SIG", "ENCFF2PEA", "2", "lab A", "cell A"],
        ["ENCSR000AAB", "ENCFF3FDR", "ENCFF3SIG", "ENCFF3PEA", "1", "lab B", "cell B"],
    ]
