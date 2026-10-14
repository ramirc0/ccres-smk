"""Build the integration test inputs: real data cut to one window per genome.

Run once from the repository root, after the frozen inputs exist:

    pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/make_fixtures.py

The lists are synthetic. Their rows point to real ENCODE files, and their
biosample names are chosen so that every branch of step 7 runs: toy_tf (DNase
and TF), toy_dnase (DNase, partial marks) and toy_nodnase (no DNase, an NA
H3K4me3 file). A toy metadata snapshot gives each name one biosample, so the
biosample keys group the rows by name. Every listed file is released in one
default analysis, so the refreshed lists keep their rows. The ENCODE mirror at
/data/projects/encode/data must exist.
"""

import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import blake3
import polars as pl

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
import freeze_input
import tf_list

MIRROR = Path("/data/projects/encode/data")
PREVIOUS = Path(
    "/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/Previous-Versions"
)
WINDOWS = {
    "hg38": ("chr21", 33_000_000, 33_500_000),
    "mm10": ("chr19", 5_000_000, 5_500_000),
}
TF_EXPERIMENTS = 6
MARK_LISTS = [f"{m}-List.txt" for m in ("DNase", "H3K4me3", "H3K27ac", "CTCF", "ATAC")]
# Every kind of freeze_input but copy, with the frozen inputs each reads.
KINDS = {
    "Hotspot-List.txt": ("hotspot-list", ["DNase-List.txt"]),
    **dict.fromkeys(MARK_LISTS, ("refresh", [])),
    "TF-List-All.txt": ("tf-list", []),
    "Biosample-Keys.txt": ("biosample-keys", [*MARK_LISTS[:4], "TF-List-All.txt"]),
    "cCREs-Previous.bed": ("previous-ccres", []),
}
ANALYSIS = "/analyses/ENCAN000TOY/"
# Output type and format of the files a hotspot row names, in column order.
HOTSPOT_FILES = list(freeze_input.HOTSPOT_FILES.items())


def rows(path):
    """Read a tab-separated file.

    Parameters
    ----------
    path : pathlib.Path
        File to read.

    Returns
    -------
    list of list of str
        Fields per line.
    """
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f]


def write(path, table):
    """Write rows as tab-separated lines, creating the directory.

    Parameters
    ----------
    path : pathlib.Path
        File to write.
    table : list of list of str
        Rows.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.writelines("\t".join(r) + "\n" for r in table)


def in_window(fields, window):
    """Tell whether a BED row lies inside the window.

    Parameters
    ----------
    fields : list of str
        BED fields.
    window : tuple of (str, int, int)
        Chromosome, start and end.

    Returns
    -------
    bool
        True if the row starts and ends within the window.
    """
    chrom, start, end = window
    return fields[0] == chrom and int(fields[1]) >= start and int(fields[2]) <= end


def cut_bed(source, target, window, compress):
    """Copy the rows of a BED file that lie in the window.

    Parameters
    ----------
    source : pathlib.Path
        BED file, gzipped if its name ends in .gz.
    target : pathlib.Path
        BED file to write.
    window : tuple of (str, int, int)
        Chromosome, start and end.
    compress : bool
        Gzip the output.
    """
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, "rt") as f:
        text = "".join(line for line in f if in_window(line.split("\t"), window))
    target.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 keeps the gzip bytes the same from run to run.
    data = text.encode()
    target.write_bytes(gzip.compress(data, mtime=0) if compress else data)


def cut_bigwig(source, target, window, sizes):
    """Cut a bigWig to the window with the UCSC tools.

    Parameters
    ----------
    source, target : pathlib.Path
        Input and output bigWig.
    window : tuple of (str, int, int)
        Chromosome, start and end.
    sizes : pathlib.Path
        Chromosome sizes.
    """
    chrom, start, end = window
    target.parent.mkdir(parents=True, exist_ok=True)
    graph = target.with_suffix(".bedGraph")
    subprocess.run(
        [
            "bigWigToBedGraph",
            f"-chrom={chrom}",
            f"-start={start}",
            f"-end={end}",
            source,
            graph,
        ],
        check=True,
    )
    subprocess.run(["bedGraphToBigWig", graph, sizes, target], check=True)
    graph.unlink()


def build(genome):
    """Write the fixture lists and files of one genome.

    Parameters
    ----------
    genome : str
        hg38 or mm10.

    Returns
    -------
    tuple of (dict, list of tuple of str, list of list of str)
        Frozen input name to fixture source path, relative to the test dir,
        the experiment, file and biosample name of every list row, and the
        hotspot list rows.
    """
    window = WINDOWS[genome]
    frozen = ROOT / "data" / "inputs" / genome
    res = HERE / "resources" / genome
    encode = HERE / "resources" / "encode"
    sizes = res / "chromInfo.txt"
    sizes.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(frozen / "chromInfo.txt", sizes)

    def encode_copy(exp, acc, ext):
        source, target = MIRROR / exp / f"{acc}.{ext}", encode / exp / f"{acc}.{ext}"
        if ext == "bigWig":
            cut_bigwig(source, target, window, sizes)
        else:
            cut_bed(source, target, window, compress=True)

    hotspots = rows(frozen / "Hotspot-List.txt")[:3]
    for r in hotspots:
        encode_copy(r[0], r[1], "bed.gz")
        encode_copy(r[0], r[2], "bigWig")
    write(res / "Hotspot-List.txt", hotspots)

    # A mark list holds one file per experiment, as in the real lists.
    first = {}
    for r in hotspots:
        first.setdefault(r[0], r)
    marks = {}
    marks["DNase"] = [
        [r[0], r[2], name] for r, name in zip(first.values(), ("toy_tf", "toy_dnase"))
    ]
    for mark, names in (
        ("H3K4me3", ("toy_tf",)),
        ("H3K27ac", ("toy_tf", "toy_dnase")),
        ("CTCF", ("toy_tf", "toy_nodnase")),
        ("ATAC", ("toy_atac",)),
    ):
        real = rows(frozen / f"{mark}-List.txt")[: len(names)]
        marks[mark] = [[r[0], r[1], name] for r, name in zip(real, names)]
    marks["H3K4me3"].append(["ENCSR999NAX", "NA", "toy_nodnase"])
    for mark, table in marks.items():
        for exp, acc, _ in table:
            if acc != "NA":
                encode_copy(exp, acc, "bigWig")
        write(res / f"{mark}-List.txt", table)

    tf_all = rows(frozen / "TF-List-All.txt")
    passing = [r for r in tf_all if tf_list.passes(r)]
    counts = []
    for r in passing[:200]:
        with gzip.open(MIRROR / r[0] / f"{r[1]}.bed.gz", "rt") as f:
            counts.append((sum(in_window(line.split("\t"), window) for line in f), r))
    chosen = [r for _, r in sorted(counts, key=lambda c: -c[0])[:TF_EXPERIMENTS]]
    for r in chosen:
        encode_copy(r[0], r[1], "bed.gz")
    write(res / "TF-List-All.txt", chosen + [r for r in tf_all if r[1] == "NA"][:1])
    named = [row for table in marks.values() for row in table]
    named += [
        (r[0], r[1], "toy_tf" if r in chosen[:2] else "toy_tf_other") for r in chosen
    ]

    # cDHSs cover only the first half, so that some rDHSs fall in the none class.
    chrom, start, end = window
    cut_bed(
        frozen / "cDHSs.bed",
        res / "cDHSs.bed",
        (chrom, start, (start + end) // 2),
        False,
    )
    for name in ("TSS-Basic.bed", "TSS-Basic-4K.bed"):
        cut_bed(frozen / name, res / name, window, compress=False)
    cut_bed(
        PREVIOUS / f"rDHS.{genome}.bed.gz", res / "rDHS-Previous.bed.gz", window, True
    )
    cut_bed(PREVIOUS / f"cCREs.{genome}.bed.gz", res / "cCREs.bed.gz", window, True)
    shutil.copy(frozen / "cres.as", res / "cres.as")

    rel = f"resources/{genome}"
    sources = {
        name: f"{rel}/{name}"
        for name in (
            "Hotspot-List.txt",
            *(f"{m}-List.txt" for m in marks),
            "TF-List-All.txt",
            "cDHSs.bed",
            "rDHS-Previous.bed.gz",
            "TSS-Basic.bed",
            "TSS-Basic-4K.bed",
            "chromInfo.txt",
            "cres.as",
        )
    }
    sources["cCREs-Previous.bed"] = f"{rel}/cCREs.bed.gz"
    return sources, named, hotspots


def write_metadata(built, meta):
    """Write a metadata snapshot that gives each biosample name one biosample.

    A hotspot row's files form its replicate. Each other list file gets a
    replicate and library of its own. Every file is released in its
    experiment's default analysis. A hotspot experiment keeps its listed lab
    and biosample summary.

    Parameters
    ----------
    built : dict
        Output of `build` per genome.
    meta : pathlib.Path
        Directory to write the Parquet files in.
    """
    meta.mkdir(parents=True, exist_ok=True)
    assembly = assemblies()
    files, experiments, listed = {}, {}, set()
    for genome, (_, named, hotspots) in built.items():
        for r in hotspots:
            experiments[r[0]] = (r[5], r[6])
            for acc, (output, form) in zip(r[1:4], HOTSPOT_FILES):
                files[acc] = (r[0], r[4], output, form, assembly[genome])
        for exp, acc, name in sorted({tuple(row) for row in named}):
            experiments.setdefault(exp, ("lab", name))
            if acc == "NA":
                continue
            listed.add((exp, acc, name))
            if acc not in files:
                rep = str(1 + sum(f[0] == exp for f in files.values()))
                files[acc] = (exp, rep, None, None, assembly[genome])
    listed = sorted(listed)
    names = sorted({name for _, _, name in listed})
    string = pl.String
    pl.DataFrame(
        {
            "id": [f"/biosamples/{n}/" for n in names],
            "biosample_ontology": [
                f"/biosample-types/toy_TOY_{i}/" for i in range(len(names))
            ],
            "donor": "/human-donors/ENCDO000TOY/",
            **dict.fromkeys(
                [
                    "age",
                    "age_units",
                    "life_stage",
                    "sex",
                    "treatments",
                    "genetic_modifications",
                    "originated_from",
                    "origin_batch",
                    "post_differentiation_time_units",
                    "disease_ontology_ids",
                    "cell_cycle_phase",
                ],
                pl.Series([None] * len(names), dtype=string),
            ),
            "post_differentiation_time": pl.Series([None] * len(names), dtype=pl.Int64),
        }
    ).write_parquet(meta / "biosample.parquet")
    library = [f"/libraries/{exp}-{files[acc][1]}/" for exp, acc, _ in listed]
    pl.DataFrame(
        {"id": library, "biosample": [f"/biosamples/{name}/" for _, _, name in listed]}
    ).write_parquet(meta / "library.parquet")
    pl.DataFrame(
        {
            "experiment": [f"/experiments/{exp}/" for exp, _, _ in listed],
            "biological_replicate": [int(files[acc][1]) for _, acc, _ in listed],
            "library": library,
            "status": "released",
        }
    ).write_parquet(meta / "replicate.parquet")
    pl.DataFrame(
        [(exp, "2020-01-01", ANALYSIS, *e) for exp, e in sorted(experiments.items())],
        schema=[
            "accession",
            "date_released",
            "default_analysis",
            "lab",
            "biosample_summary",
        ],
        orient="row",
    ).write_parquet(meta / "experiment.parquet")
    pl.DataFrame(
        [
            (acc, f"/experiments/{exp}/", rep, output, form, genome, ANALYSIS)
            for acc, (exp, rep, output, form, genome) in sorted(files.items())
        ],
        schema={
            "accession": string,
            "dataset": string,
            "biological_replicates": string,
            "output_type": string,
            "file_format": string,
            "genome_assembly": string,
            "analyses": string,
        },
        orient="row",
    ).with_columns(status=pl.lit("released")).write_parquet(meta / "file.parquet")
    treatment = ["id", "biosample_treatment_term_id", "biosample_treatment"]
    for name, schema in (
        (
            "treatment",
            {
                **dict.fromkeys(treatment, string),
                "amount": pl.Float64,
                "amount_units": string,
                "duration": pl.Float64,
                "duration_units": string,
            },
        ),
        ("genetic_modification", dict.fromkeys(["id", "accession", "purpose"], string)),
        ("mouse_donor", dict.fromkeys(["id", "strain_background"], string)),
        (
            "chip_peak_enrichment_quality_metric",
            {"files": string, "status": string, "frip_score_of_peaks": pl.Float64},
        ),
    ):
        pl.DataFrame(schema=schema).write_parquet(meta / f"{name}.parquet")


def template_head():
    """Return the config template up to its inputs.

    Returns
    -------
    str
        Run settings and per-genome constants.
    """
    config = (ROOT / "config" / "config.yaml.template").read_text()
    return config[: config.index("\ninputs:\n") + 1]


def assemblies():
    """Read each genome's ENCODE assembly name from the config template.

    Returns
    -------
    dict of str to str
        Assembly per genome.
    """
    return dict(
        re.findall(r"^  (\w+):\n    assembly: (\S+)", template_head(), re.MULTILINE)
    )


def entries(sources):
    """List a genome's frozen inputs, each with how the test config builds it.

    Parameters
    ----------
    sources : dict of str to str
        Frozen input name to fixture source path.

    Returns
    -------
    dict of str to tuple
        Kind, source path or None, and the frozen inputs it reads.
    """
    specs = {}
    for name in [*sources, "Biosample-Keys.txt"]:
        kind, frozen = KINDS.get(name, ("copy", []))
        specs[name] = (kind, sources.get(name), frozen)
    return specs


def main():
    """Build both genomes' fixtures, then the test config and checksums."""
    built = {genome: build(genome) for genome in WINDOWS}
    meta = HERE / "resources" / "metadata"
    write_metadata(built, meta)
    head = template_head()
    head = head.replace(
        "encode_root: /data/projects/encode/data", "encode_root: resources/encode"
    )
    head = re.sub(r"encode_metadata: \S+", "encode_metadata: resources/metadata", head)
    lines = [head, "inputs:\n"]
    sums = []
    assembly = assemblies()
    (HERE / "data").mkdir(exist_ok=True)
    for genome, (sources, _, _) in built.items():
        lines.append(f"  {genome}:\n")
        specs = entries(sources)
        for name, (kind, src, frozen) in specs.items():
            if kind == "copy":
                lines.append(f"    {name}: {src}\n")
                continue
            lines.append(f"    {name}:\n      kind: {kind}\n")
            if src:
                lines.append(f"      source: {src}\n")
            if frozen:
                lines.append(f"      frozen: [{', '.join(frozen)}]\n")
        digests = {}
        with tempfile.TemporaryDirectory(dir=HERE / "data") as tmp:
            # Inputs that read frozen ones come after them.
            for name in sorted(specs, key=lambda n: bool(specs[n][2])):
                kind, src, frozen = specs[name]
                out = Path(tmp, name)
                freeze_input.build(
                    kind,
                    [*([HERE / src] if src else []), *(Path(tmp, f) for f in frozen)],
                    meta,
                    assembly[genome],
                    out,
                )
                digests[name] = blake3.blake3(out.read_bytes()).hexdigest()
        sums += [f"{digests[name]}  {genome}/{name}\n" for name in specs]
    (HERE / "config").mkdir(exist_ok=True)
    (HERE / "config" / "config.yaml").write_text("".join(lines))
    (HERE / "config" / "inputs.b3").write_text("".join(sums))


if __name__ == "__main__":
    main()
