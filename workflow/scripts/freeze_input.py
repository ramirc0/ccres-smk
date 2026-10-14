"""Freeze one pipeline input: copy or derive it, then check its BLAKE3 sum.

The sources sit in another user's tree and keep changing. A frozen copy whose
checksum is not recorded in the checksum file fails the job.
"""

import argparse
import gzip
import shutil
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
    p.add_argument(
        "--sources",
        nargs="+",
        required=True,
        help="Source file, then the frozen inputs the kind reads.",
    )
    p.add_argument(
        "--kind",
        default="copy",
        choices=[
            "copy",
            "previous-ccres",
            "refresh",
            "tf-list",
            "hotspot-list",
            "biosample-keys",
        ],
        help="How to build the frozen file from its sources.",
    )
    p.add_argument(
        "--metadata", type=Path, help="ENCODE metadata snapshot, for derived lists."
    )
    p.add_argument("--assembly", help="Genome assembly, as ENCODE names it.")
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
    meta : pathlib.Path
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
    experiment_keys(pairs, meta).write_csv(
        out, separator="\t", include_header=False, quote_style="never"
    )


FILE_COLUMNS = [
    "accession",
    "dataset",
    "status",
    "output_type",
    "file_format",
    "biological_replicates",
    "genome_assembly",
    "analyses",
]
# The DNase pipeline's per-replicate outputs that a hotspot list row names.
HOTSPOT_FILES = {
    "FDR cut rate": "bed",
    "read-depth normalized signal": "bigWig",
    "peaks": "starch",
}
# No analysis holds both names. Both pair with conservative IDR peaks.
RENAMED = {"optimal IDR thresholded peaks": "IDR thresholded peaks"}


def file_table(meta):
    """Read the file report of a metadata snapshot, each file with its experiment.

    Parameters
    ----------
    meta : pathlib.Path
        Directory of the ENCODE metadata snapshot.

    Returns
    -------
    polars.DataFrame
        The `FILE_COLUMNS` and exp, the accession of the file's dataset now.
    """
    return pl.read_parquet(meta / "file.parquet", columns=FILE_COLUMNS).with_columns(
        exp=pl.col("dataset").str.extract(r"^/experiments/(\w+)/$")
    )


def default_files(files, meta, assembly):
    """Keep the released files of each experiment's default analysis.

    Parameters
    ----------
    files : polars.DataFrame
        Output of `file_table`.
    meta : pathlib.Path
        Directory of the ENCODE metadata snapshot.
    assembly : str
        Genome assembly, as ENCODE names it (GRCh38, mm10).

    Returns
    -------
    polars.DataFrame
        The released files of the assembly in their experiment's default
        analysis.
    """
    analysis = pl.read_parquet(
        meta / "experiment.parquet", columns=["accession", "default_analysis"]
    )
    return (
        files.filter(
            (pl.col("status") == "released") & (pl.col("genome_assembly") == assembly)
        )
        .join(analysis, left_on="exp", right_on="accession")
        .filter(
            pl.col("analyses").str.contains(pl.col("default_analysis"), literal=True)
        )
    )


def refresh(rows, meta, assembly):
    """Drop list rows whose file was withdrawn and swap archived files.

    A kept row takes its file's experiment now. An archived file gives way to
    the released file of its experiment's default analysis with the same
    format and output type whose replicates are the largest subset of its
    own. ENCODE4 calls the optimal IDR thresholded peaks of older pipelines
    `IDR thresholded peaks`. A file neither released nor archived, such as a
    revoked one, drops its row. A file of NA stays as it is.

    Parameters
    ----------
    rows : list of list of str
        List rows: experiment, file, then columns that stay.
    meta : pathlib.Path
        Directory of the ENCODE metadata snapshot.
    assembly : str
        Genome assembly, as ENCODE names it.

    Returns
    -------
    list of list of str
        The refreshed rows, in order.

    Raises
    ------
    ValueError
        If an archived file has no single replacement.
    """
    files = file_table(meta).with_columns(
        pl.col("output_type").replace(RENAMED),
        replicates=pl.col("biological_replicates").str.split(","),
    )
    listed = files.filter(pl.col("accession").is_in([row[1] for row in rows]))
    same = ["exp", "output_type", "file_format"]
    archived = listed.filter(pl.col("status") == "archived").select(
        *same, old="accession", old_replicates="replicates"
    )
    size = pl.col("replicates").list.len()
    found = (
        default_files(files, meta, assembly)
        .join(archived, on=same)
        .filter(
            pl.col("replicates").list.set_difference("old_replicates").list.len() == 0
        )
        .group_by("old")
        .agg(pl.col("accession").filter(size == size.max()))
    )
    swaps = {old: new[0] for old, new in found.iter_rows() if len(new) == 1}
    if missing := set(archived["old"]) - set(swaps):
        raise ValueError(f"No single replacement for the archived files {missing}.")
    records = files.filter(
        pl.col("accession").is_in([*listed["accession"], *swaps.values()])
    )
    status = dict(records.select("accession", "status").iter_rows())
    experiment = dict(records.select("accession", "exp").iter_rows())
    out = []
    for exp, acc, *rest in rows:
        acc = swaps.get(acc, acc)
        if acc == "NA":
            out.append([exp, acc, *rest])
        elif status[acc] == "released":
            out.append([experiment[acc], acc, *rest])
    return out


def encode_frip(rows, meta):
    """Set each TF list row's FRiP to ENCODE's FRiP of its peak file.

    The listed FRiPs are ENCODE's peak enrichment QC, where a file has one. A
    row whose file has none keeps its FRiP.

    Parameters
    ----------
    rows : list of list of str
        TF list rows: experiment, peak file, biosample, target, FRiP, then
        other columns.
    meta : pathlib.Path
        Directory of the ENCODE metadata snapshot.

    Returns
    -------
    list of list of str
        The rows, in order.
    """
    frip = dict(
        pl.read_parquet(
            meta / "chip_peak_enrichment_quality_metric.parquet",
            columns=["files", "status", "frip_score_of_peaks"],
        )
        .filter(pl.col("status") == "released")
        .select(pl.col("files").str.extract(r"^/files/(\w+)/$"), "frip_score_of_peaks")
        .drop_nulls()
        .iter_rows()
    )
    return [[*row[:4], str(frip.get(row[1], row[4])), *row[5:]] for row in rows]


def hotspot_rows(hotspots, dnase, meta, assembly):
    """Rebuild the hotspot list from the default analyses of its experiments.

    A row is one biological replicate: the FDR cut rate, normalized signal
    and Hotspot2 peaks it has in the default analysis of the experiment that
    holds them now, then the replicate, the lab and the biosample summary.
    The listed rows come first, in order, each rebuilt from the experiment and
    replicate of its FDR file. A row whose FDR file was withdrawn is dropped.
    The other replicates of the listed experiments follow, then those of the
    DNase list's experiments that the hotspot list lacks.

    Parameters
    ----------
    hotspots : list of list of str
        Hotspot list rows: experiment, FDR file, then other columns.
    dnase : list of list of str
        DNase list rows: experiment, signal file, then other columns.
    meta : pathlib.Path
        Directory of the ENCODE metadata snapshot.
    assembly : str
        Genome assembly, as ENCODE names it.

    Returns
    -------
    list of list of str
        Experiment, FDR file, signal file, peaks file, replicate, lab and
        biosample summary of each row.

    Raises
    ------
    ValueError
        If a listed replicate has no row, or two files of one kind.
    """
    files = file_table(meta)
    listed = {
        f["accession"]: f
        for f in files.filter(
            pl.col("accession").is_in([r[1] for r in hotspots + dnase])
        ).iter_rows(named=True)
        if f["status"] in ("released", "archived")
    }
    keys = [
        (listed[r[1]]["exp"], listed[r[1]]["biological_replicates"])
        for r in hotspots
        if r[1] in listed
    ]
    dnase_exps = [listed[r[1]]["exp"] for r in dnase if r[1] in listed]
    order = list(dict.fromkeys([e for e, _ in keys] + dnase_exps))
    wanted = [{"output_type": o, "file_format": f} for o, f in HOTSPOT_FILES.items()]
    found = default_files(files, meta, assembly).filter(
        pl.col("exp").is_in(order)
        & pl.struct("output_type", "file_format").is_in(wanted)
        & ~pl.col("biological_replicates").str.contains(",")
    )
    kind = ["exp", "biological_replicates", "output_type"]
    if found.select(kind).is_duplicated().any():
        raise ValueError("A replicate has two files of one kind in its analysis.")
    experiments = pl.read_parquet(
        meta / "experiment.parquet", columns=["accession", "lab", "biosample_summary"]
    )
    table = (
        found.group_by("exp", "biological_replicates")
        .agg(
            pl.col("accession")
            .filter(pl.col("output_type") == output)
            .first()
            .alias(output)
            for output in HOTSPOT_FILES
        )
        .drop_nulls()
        .join(experiments, left_on="exp", right_on="accession")
    )
    rows = {
        (exp, rep): [exp, fdr, signal, peaks, rep, lab, biosample]
        for exp, rep, fdr, signal, peaks, lab, biosample in table.select(
            "exp", "biological_replicates", *HOTSPOT_FILES, "lab", "biosample_summary"
        ).iter_rows()
    }
    if missing := [k for k in keys if k not in rows]:
        raise ValueError(f"No hotspot files for the listed replicates {missing}.")
    for exp in order:
        keys += sorted((k for k in rows if k[0] == exp), key=lambda k: int(k[1]))
    return [rows[k] for k in dict.fromkeys(keys)]


def read_rows(path):
    """Read a tab-separated list.

    Parameters
    ----------
    path : str
        List file.

    Returns
    -------
    list of list of str
        Fields of each line.
    """
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f]


def write_rows(rows, path):
    """Write rows as tab-separated lines.

    Parameters
    ----------
    rows : list of list of str
        Fields of each line.
    path : str
        List file to write.
    """
    with open(path, "w") as f:
        f.writelines("\t".join(row) + "\n" for row in rows)


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


def build(kind, sources, meta, assembly, out):
    """Build a frozen input from its sources.

    Parameters
    ----------
    kind : str
        How: copy, previous-ccres, refresh, tf-list, hotspot-list or
        biosample-keys.
    sources : list of str
        Source file, then the frozen inputs the kind reads.
    meta : pathlib.Path or None
        Directory of the ENCODE metadata snapshot, for derived lists.
    assembly : str or None
        Genome assembly, as ENCODE names it, for refreshed lists.
    out : str
        File to write.
    """
    if kind == "copy":
        (src,) = sources
        shutil.copyfile(src, out)
    elif kind == "previous-ccres":
        (src,) = sources
        previous_ccres(src, out)
    elif kind in ("refresh", "tf-list"):
        (src,) = sources
        rows = refresh(read_rows(src), meta, assembly)
        if kind == "tf-list":
            rows = encode_frip(rows, meta)
        write_rows(rows, out)
    elif kind == "hotspot-list":
        hotspots, dnase = map(read_rows, sources)
        write_rows(hotspot_rows(hotspots, dnase, meta, assembly), out)
    else:
        write_biosample_keys(sources, meta, out)


def main(argv=None):
    """Build the frozen file and fail unless its checksum is the recorded one.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    build(args.kind, args.sources, args.metadata, args.assembly, args.output)

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
