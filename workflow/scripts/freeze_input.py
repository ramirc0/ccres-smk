"""Freeze one pipeline input: copy or derive it, then check its BLAKE3 sum.

The sources sit in another user's tree and keep changing. A frozen copy whose
checksum is not recorded in the checksum file fails the job.
"""

import argparse
import gzip
import shutil
from contextlib import ExitStack
from pathlib import Path

import blake3
import polars as pl


def build_parser():
    """Build the freeze_input.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the freeze_input.py flags.
    """
    p = argparse.ArgumentParser(description="Freeze one pipeline input.")
    p.add_argument("--sources", nargs="+", required=True, help="Source file(s).")
    p.add_argument(
        "--kind",
        default="copy",
        choices=["copy", "previous-ccres", "cell-type-tf", "biosample-keys"],
        help="How to build the frozen file from its sources.",
    )
    p.add_argument(
        "--checksums", required=True, help="b3sum-format file of expected sums."
    )
    p.add_argument(
        "--key", required=True, help="Name of the input in the checksum file."
    )
    p.add_argument("-o", "--output", required=True, help="Frozen file to write.")
    return p


def previous_ccres(src, out):
    """Write the V1-V3 rows of the previous-versions cCRE file.

    The source now also holds the V4 release, built by the run this pipeline
    reproduces. Its accessions must not seed the V4 accessioning.

    Parameters
    ----------
    src : str
        Gzipped cCRE BED. Column 6 is the version that first reported the cCRE.
    out : str
        BED to write, uncompressed.
    """
    with gzip.open(src, "rt") as f, open(out, "w") as o:
        o.writelines(line for line in f if line.rstrip("\n").split("\t")[5] != "V4")


def cell_type_tf(tf_columns, master, out):
    """Write the TF experiments of the biosamples flagged TF=yes in May 2022.

    The step 7 TF list kept growing after the release. Its rows for the
    biosamples the release flagged keep the release's TF calls reproducible.

    Parameters
    ----------
    tf_columns : list of str
        TF experiment list, as files whose lines are joined side by side with
        tabs. The last column is the biosample.
    master : str
        Master cell list of the release. Column 9 is the biosample, the last
        column `yes` or `no`.
    out : str
        Rows of the TF list whose biosample is flagged `yes`, in file order.
    """
    with open(master) as f:
        flagged = {
            fields[8]
            for fields in (line.rstrip("\n").split("\t") for line in f)
            if fields[-1] == "yes"
        }
    with ExitStack() as stack, open(out, "w") as o:
        files = [stack.enter_context(open(path)) for path in tf_columns]
        for parts in zip(*files, strict=True):
            line = "\t".join(part.rstrip("\n") for part in parts)
            if line.split("\t")[-1] in flagged:
                o.write(f"{line}\n")


def _ids(column):
    """Split a column of comma-separated @ids into lists.

    Parameters
    ----------
    column : str
        Column name.

    Returns
    -------
    polars.Expr
        List of the @ids, empty for a null.
    """
    ids = pl.col(column).fill_null("").str.split(",")
    return ids.list.eval(pl.element().filter(pl.element() != ""))


def _field(name, value):
    """Return `name=value`, or null when the value is null or empty.

    Parameters
    ----------
    name : str
        Field name.
    value : polars.Expr
        Field value, a string.

    Returns
    -------
    polars.Expr
        The field text.
    """
    return pl.when(value.is_not_null() & (value != "")).then(
        pl.format(f"{name}={{}}", value)
    )


def _key(fields):
    """Join the present fields of a key with `|`.

    Parameters
    ----------
    fields : dict of str to polars.Expr
        Field values by name, in key order.

    Returns
    -------
    polars.Expr
        The key text.
    """
    return pl.concat_str(
        [_field(name, value) for name, value in fields.items()],
        separator="|",
        ignore_nulls=True,
    )


def biosample_keys(meta):
    """Key each biosample by the attributes that tell two samples apart.

    The key joins, in order: ontology term, donor (a mouse's strain
    background), treatments (term, amount and duration, sorted), age, life
    stage, sex when male or female, post-differentiation time, genetic
    modifications, cell-cycle phase, diseases, and the term the sample
    originated from when it differs. A missing post-differentiation time comes
    from the samples of the same origin batch and term, when they agree on one.

    Parameters
    ----------
    meta : pathlib.Path
        Directory of the ENCODE metadata snapshot, one Parquet file per type.

    Returns
    -------
    polars.DataFrame
        Biosample @id, its key, and its key without tagging modifications.
    """
    bio = pl.read_parquet(
        meta / "biosample.parquet",
        columns=[
            "id",
            "biosample_ontology",
            "donor",
            "age",
            "age_units",
            "life_stage",
            "sex",
            "treatments",
            "genetic_modifications",
            "originated_from",
            "origin_batch",
            "post_differentiation_time",
            "post_differentiation_time_units",
            "disease_ontology_ids",
            "cell_cycle_phase",
        ],
    ).with_columns(
        term=pl.col("biosample_ontology")
        .str.strip_prefix("/biosample-types/")
        .str.strip_suffix("/"),
        differentiation=pl.when(pl.col("post_differentiation_time").is_not_null()).then(
            pl.concat_str(
                "post_differentiation_time",
                "post_differentiation_time_units",
                separator=" ",
                ignore_nulls=True,
            )
        ),
    )
    batch = (
        bio.drop_nulls(["origin_batch", "differentiation"])
        .group_by("origin_batch", "term")
        .agg(pl.col("differentiation").unique())
        .filter(pl.col("differentiation").list.len() == 1)
        .select("origin_batch", "term", batch=pl.col("differentiation").list.first())
    )
    treatment = pl.read_parquet(meta / "treatment.parquet").select(
        "id",
        text=pl.concat_str(
            pl.coalesce("biosample_treatment_term_id", "biosample_treatment"),
            pl.format("{} {}", "amount", "amount_units"),
            pl.format("for {} {}", "duration", "duration_units"),
            separator=" ",
            ignore_nulls=True,
        ),
    )
    treatments = (
        bio.select("id", t=_ids("treatments"))
        .explode("t")
        .join(treatment, left_on="t", right_on="id")
        .group_by("id")
        .agg(treatment_text=pl.col("text").sort().str.join(","))
    )
    modification = pl.read_parquet(
        meta / "genetic_modification.parquet", columns=["id", "accession", "purpose"]
    )
    modifications = (
        bio.select("id", m=_ids("genetic_modifications"))
        .explode("m")
        .join(modification, left_on="m", right_on="id")
        .group_by("id")
        .agg(
            modifications=pl.col("accession").sort().str.join(","),
            untagged=pl.col("accession")
            .filter(pl.col("purpose").ne_missing("tagging"))
            .sort()
            .str.join(","),
        )
    )
    strains = pl.read_parquet(
        meta / "mouse_donor.parquet", columns=["id", "strain_background"]
    )
    bio = (
        bio.join(batch, on=["origin_batch", "term"], how="left")
        .join(treatments, on="id", how="left")
        .join(modifications, on="id", how="left")
        .join(strains, left_on="donor", right_on="id", how="left")
        .join(
            bio.select(originated_from="id", origin="term"),
            on="originated_from",
            how="left",
        )
    )
    donor = (
        pl.when(pl.col("donor").str.starts_with("/mouse-donors/"))
        .then(pl.col("strain_background"))
        .otherwise(pl.col("donor").str.extract(r"/([^/]+)/$"))
    )
    fields = {
        "term": pl.col("term"),
        "donor": donor,
        "treatments": pl.col("treatment_text"),
        "age": pl.concat_str("age", "age_units", separator=" ", ignore_nulls=True),
        "life_stage": pl.col("life_stage"),
        "sex": pl.when(pl.col("sex").is_in(["male", "female"])).then(pl.col("sex")),
        "differentiation": pl.coalesce("differentiation", "batch"),
        "modifications": pl.col("modifications"),
        "cell_cycle": pl.col("cell_cycle_phase"),
        "disease": _ids("disease_ontology_ids").list.sort().list.join(","),
        "origin": pl.when(pl.col("origin") != pl.col("term")).then(pl.col("origin")),
    }
    return bio.select(
        "id",
        key=_key(fields),
        tf_key=_key({**fields, "modifications": pl.col("untagged")}),
    )


def experiment_keys(pairs, meta):
    """Key each listed experiment file by the biosamples it came from.

    A file's biosamples are those of the libraries of its dataset's
    replicates, among the biological replicates it lists. A file that lists
    none, as revoked ones do, takes every replicate of its dataset.

    Parameters
    ----------
    pairs : polars.DataFrame
        Columns exp and file, as the lists name them.
    meta : pathlib.Path
        Directory of the ENCODE metadata snapshot.

    Returns
    -------
    polars.DataFrame
        Columns exp, file, key, tf_key and the experiment's release date, the
        keys being the sorted biosample keys joined by " + ", sorted by exp
        and file.

    Raises
    ------
    ValueError
        If a file has no biosample.
    """
    files = pl.read_parquet(
        meta / "file.parquet", columns=["accession", "dataset", "biological_replicates"]
    )
    replicates = pl.read_parquet(
        meta / "replicate.parquet",
        columns=["experiment", "biological_replicate", "library", "status"],
    ).filter(pl.col("status") != "deleted")
    libraries = pl.read_parquet(meta / "library.parquet", columns=["id", "biosample"])
    released = pl.read_parquet(
        meta / "experiment.parquet", columns=["accession", "date_released"]
    )
    numbers = pl.col("biological_replicates").str.split(",").cast(pl.List(pl.Int64))
    joined = (
        pairs.join(files, left_on="file", right_on="accession", how="left")
        .with_columns(numbers=numbers)
        .join(replicates, left_on="dataset", right_on="experiment", how="left")
        .filter(
            pl.col("numbers").is_null()
            | pl.col("numbers").list.contains(pl.col("biological_replicate"))
        )
        .join(libraries, left_on="library", right_on="id", how="left")
        .join(biosample_keys(meta), left_on="biosample", right_on="id", how="left")
    )
    keyed = (
        joined.drop_nulls("key")
        .group_by("exp", "file")
        .agg(pl.col(k).unique().sort().str.join(" + ") for k in ("key", "tf_key"))
    )
    missing = pairs.join(keyed, on=["exp", "file"], how="anti")
    if missing.height:
        raise ValueError(f"No biosample for {missing.rows()}.")
    return keyed.join(released, left_on="exp", right_on="accession", how="left").sort(
        "exp", "file"
    )


def write_biosample_keys(lists, meta, out):
    """Write the biosample keys of the experiment files that lists name.

    Parameters
    ----------
    lists : list of str
        Lists whose first two columns are experiment and file. A file of NA
        is skipped.
    meta : str
        Directory of the ENCODE metadata snapshot.
    out : str
        Table to write: experiment, file, key, key without tagging
        modifications, and release date.
    """
    pairs = set()
    for path in lists:
        with open(path) as f:
            pairs |= {tuple(line.split("\t")[:2]) for line in f}
    pairs = pl.DataFrame(
        sorted(p for p in pairs if p[1] != "NA"), schema=["exp", "file"], orient="row"
    )
    experiment_keys(pairs, Path(meta)).write_csv(
        out, separator="\t", include_header=False, quote_style="never"
    )


def read_checksums(path):
    """Read a b3sum-format checksum file.

    Parameters
    ----------
    path : str
        Lines of `<hex digest>  <key>`.

    Returns
    -------
    dict of str to str
        Digest per key.
    """
    with open(path) as f:
        return {key: digest for digest, key in (line.split() for line in f)}


def main(argv=None):
    """Build the frozen file and fail unless its checksum is the recorded one.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    if args.kind == "copy":
        (src,) = args.sources
        shutil.copyfile(src, args.output)
    elif args.kind == "previous-ccres":
        (src,) = args.sources
        previous_ccres(src, args.output)
    elif args.kind == "biosample-keys":
        *lists, meta = args.sources
        write_biosample_keys(lists, meta, args.output)
    else:
        *tf_columns, master = args.sources
        cell_type_tf(tf_columns, master, args.output)

    digest = blake3.blake3(max_threads=blake3.blake3.AUTO).update_mmap(args.output)
    digest = digest.hexdigest()
    expected = read_checksums(args.checksums).get(args.key)
    if expected is None:
        raise SystemExit(f"{args.key}: no checksum in {args.checksums}. Got {digest}.")
    if digest != expected:
        raise SystemExit(f"{args.key}: checksum {digest}, expected {expected}.")
    print(f"{args.key}: {digest} OK")


if __name__ == "__main__":
    main()
