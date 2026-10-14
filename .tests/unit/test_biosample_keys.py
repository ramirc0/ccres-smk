"""Tests of the biosample keys that freeze_input derives from ENCODE metadata."""

import freeze_input
import polars as pl

BIOSAMPLE = {
    "id": None,
    "biosample_ontology": "/biosample-types/cell_line_EFO_0002067/",
    "donor": "/human-donors/ENCDO000AAD/",
    "age": "53",
    "age_units": "year",
    "life_stage": "adult",
    "sex": "female",
    "treatments": None,
    "genetic_modifications": None,
    "originated_from": None,
    "origin_batch": None,
    "post_differentiation_time": None,
    "post_differentiation_time_units": None,
    "disease_ontology_ids": None,
    "cell_cycle_phase": None,
}


def metadata(tmp_path, biosamples, treatments=(), modifications=(), strains=()):
    """Write a metadata snapshot of the given objects.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to write the Parquet files in.
    biosamples : list of dict
        Biosample fields that differ from `BIOSAMPLE`. Each needs an id.
    treatments : sequence of tuple
        Treatment id, term id, name, amount, amount units, duration, duration
        units.
    modifications : sequence of tuple of str
        Modification id, accession, purpose.
    strains : sequence of tuple of str
        Mouse donor id, strain background.

    Returns
    -------
    pathlib.Path
        The snapshot directory.
    """
    string = pl.String
    pl.DataFrame(
        [{**BIOSAMPLE, **b} for b in biosamples],
        schema={
            **dict.fromkeys(BIOSAMPLE, string),
            "post_differentiation_time": pl.Int64,
        },
    ).write_parquet(tmp_path / "biosample.parquet")
    pl.DataFrame(
        treatments,
        schema={
            "id": string,
            "biosample_treatment_term_id": string,
            "biosample_treatment": string,
            "amount": pl.Float64,
            "amount_units": string,
            "duration": pl.Float64,
            "duration_units": string,
        },
        orient="row",
    ).write_parquet(tmp_path / "treatment.parquet")
    pl.DataFrame(
        modifications,
        schema=dict.fromkeys(["id", "accession", "purpose"], string),
        orient="row",
    ).write_parquet(tmp_path / "genetic_modification.parquet")
    pl.DataFrame(
        strains, schema=dict.fromkeys(["id", "strain_background"], string), orient="row"
    ).write_parquet(tmp_path / "mouse_donor.parquet")
    return tmp_path


def keys(tmp_path, *args, **kwargs):
    """Key the biosamples of a snapshot.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to write the snapshot in.
    *args, **kwargs
        Passed to `metadata`.

    Returns
    -------
    dict of str to tuple of (str, str)
        Key and key without tagging modifications per biosample id.
    """
    table = freeze_input.biosample_keys(metadata(tmp_path, *args, **kwargs))
    return {i: (k, t) for i, k, t in table.iter_rows()}


def test_treatment_order_does_not_matter(tmp_path):
    """Key two samples treated alike the same, whatever the treatment order."""
    found = keys(
        tmp_path,
        [
            {"id": "a", "treatments": "/t/1/,/t/2/"},
            {"id": "b", "treatments": "/t/2/,/t/1/"},
            {"id": "c", "treatments": "/t/1/"},
        ],
        treatments=[
            ("/t/1/", "CHEBI:1", "drug", 5.0, "uM", 4.0, "hour"),
            ("/t/2/", None, "beads", None, None, 24.0, "hour"),
        ],
    )
    assert found["a"] == found["b"]
    assert found["a"] != found["c"]
    assert "treatments=CHEBI:1 5.0 uM for 4.0 hour,beads for 24.0 hour" in found["a"][0]


def test_differentiation_time_splits(tmp_path):
    """Key samples differentiated for different times apart."""
    found = keys(
        tmp_path,
        [
            {"id": "a", "post_differentiation_time": 8},
            {"id": "b", "post_differentiation_time": 180},
            {"id": "c"},
        ],
    )
    assert len(set(found.values())) == 3


def test_origin_batch_fills_a_missing_time(tmp_path):
    """Take a missing time from the same-term samples of the batch."""
    batch = {"origin_batch": "/b/1/", "post_differentiation_time_units": "day"}
    found = keys(
        tmp_path,
        [
            {"id": "a", **batch},
            {"id": "b", **batch, "post_differentiation_time": 8},
            {"id": "c", **batch, "biosample_ontology": "/biosample-types/x_1/"},
        ],
    )
    assert found["a"] == found["b"]
    assert "differentiation" not in found["c"][0]


def test_mouse_donors_key_by_strain(tmp_path):
    """Key mouse samples by strain background, not by donor."""
    found = keys(
        tmp_path,
        [
            {"id": "a", "donor": "/mouse-donors/ENCDO000ABC/"},
            {"id": "b", "donor": "/mouse-donors/ENCDO000ABD/"},
            {"id": "c", "donor": "/mouse-donors/ENCDO000ABE/"},
        ],
        strains=[
            ("/mouse-donors/ENCDO000ABC/", "C57BL/6"),
            ("/mouse-donors/ENCDO000ABD/", "C57BL/6"),
            ("/mouse-donors/ENCDO000ABE/", "Castaneus"),
        ],
    )
    assert found["a"] == found["b"]
    assert found["a"] != found["c"]


def test_tagging_counts_only_in_the_full_key(tmp_path):
    """Drop tagging modifications from the TF key only."""
    found = keys(
        tmp_path,
        [
            {"id": "a", "genetic_modifications": "/gm/1/,/gm/2/"},
            {"id": "b", "genetic_modifications": "/gm/1/"},
            {"id": "c"},
        ],
        modifications=[
            ("/gm/1/", "ENCGM000AAA", "tagging"),
            ("/gm/2/", "ENCGM000AAB", None),
        ],
    )
    assert "modifications=ENCGM000AAA,ENCGM000AAB" in found["a"][0]
    assert found["a"][1].endswith("modifications=ENCGM000AAB")
    assert found["b"][0] != found["c"][0]
    assert found["b"][1] == found["c"][1]
