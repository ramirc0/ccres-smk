# Frozen inputs and ENCODE downloads.

# Derived inputs filter their sources. Every other input is a plain copy.
DERIVED = {
    "cCREs-Previous.bed": "previous-ccres",
    "Biosample-Keys.txt": "biosample-keys",
}
# Lists whose experiments Biosample-Keys.txt keys from the metadata snapshot.
KEYED_LISTS = [
    *(f"{mark}-List.txt" for mark in ("DNase", "H3K4me3", "H3K27ac", "CTCF")),
    "TF-List.All.txt",
]


def frozen_names(genome):
    """List the frozen inputs of a genome.

    Parameters
    ----------
    genome : str
        Genome key.

    Returns
    -------
    list of str
        The inputs with sources in the config, then Biosample-Keys.txt.
    """
    return [*config["sources"][genome], "Biosample-Keys.txt"]


def input_sources(w):
    """Return the source files of a frozen input.

    Sources live in another user's tree and keep changing. `ancient` stops
    their mtimes from triggering reruns. The checksum guards their content.
    Biosample-Keys.txt derives from frozen lists and the metadata snapshot.

    Parameters
    ----------
    w : snakemake.io.Wildcards
        Wildcards `genome` and `name`.

    Returns
    -------
    list of str
        Source paths, outside sources flagged ancient.
    """
    if w.name == "Biosample-Keys.txt":
        lists = [frozen(w.genome, name) for name in KEYED_LISTS]
        return [*lists, ancient(config["encode_metadata"])]
    src = config["sources"][w.genome][w.name]
    return [ancient(s) for s in ([src] if isinstance(src, str) else src)]


# Input functions read the frozen lists to enumerate jobs, so a fresh inputs_dir
# needs this target before any other.
rule freeze:
    input:
        [frozen(genome, name) for genome in GENOMES for name in frozen_names(genome)],


rule freeze_input:
    input:
        input_sources,
    output:
        f"{INPUTS}/{{genome}}/{{name}}",
    params:
        kind=lambda w: DERIVED.get(w.name, "copy"),
        key=lambda w: f"{w.genome}/{w.name}",
        checksums="config/inputs.b3",
    log:
        f"{LOGDIR}/freeze_input/{{genome}}/{{name}}.txt",
    benchmark:
        f"{BENCHDIR}/freeze_input/{{genome}}/{{name}}.tsv",
    wildcard_constraints:
        name=r"[^/]+",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/freeze_input.py \
            --sources {input:q} \
            --kind {params.kind:q} \
            --checksums {params.checksums:q} \
            --key {params.key:q} \
            --output {output:q}

        sync {output:q}
        """


def encode_files(genome):
    """List the ENCODE files a genome's jobs read.

    Parameters
    ----------
    genome : str
        Genome key.

    Returns
    -------
    list of str
        Path of each file, in the mirror or the download directory.
    """
    marks = [r for m in MARKS for r in read_rows(frozen_list(genome, f"{m}-List.txt"))]
    hotspots = read_rows(frozen_list(genome, "Hotspot-List.txt"))
    return [encode_file(*f) for f in encode_inputs(hotspots, marks, tf_rows(genome))]


# Fetches the ENCODE files the mirror lacks, so a run needs no network.
rule download:
    input:
        lambda w: [path for genome in GENOMES for path in encode_files(genome)],


# Only for files missing from the local mirror.
rule download_encode:
    output:
        f"{DOWNLOAD_DIR}/{{exp}}/{{acc}}.{{ext}}",
    log:
        f"{LOGDIR}/download_encode/{{exp}}/{{acc}}.{{ext}}.txt",
    benchmark:
        f"{BENCHDIR}/download_encode/{{exp}}/{{acc}}.{{ext}}.tsv",
    params:
        metadata=f"{config['encode_metadata']}/file.parquet",
    wildcard_constraints:
        ext=r"bed\.gz|bigWig",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/download_encode.py \
            --accession {wildcards.acc:q} \
            --metadata {params.metadata:q} \
            --output {output:q}

        sync {output:q}
        """
