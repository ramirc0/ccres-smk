# Frozen inputs and ENCODE downloads.


def freeze_spec(genome, name):
    """Return how to build a frozen input.

    Parameters
    ----------
    genome : str
        Genome key.
    name : str
        Frozen input name, a key of `config["inputs"][genome]`.

    Returns
    -------
    dict
        kind, source (an outside file, or None) and frozen (the frozen
        inputs it reads).
    """
    spec = config["inputs"][genome][name]
    if isinstance(spec, str):
        spec = {"source": spec}
    return {"kind": "copy", "source": None, "frozen": [], **spec}


def input_sources(w):
    """Return the files a frozen input is built from.

    Sources live in another user's tree and keep changing. `ancient` stops
    their mtimes from triggering reruns. The checksum guards their content.

    Parameters
    ----------
    w : snakemake.io.Wildcards
        Wildcards `genome` and `name`.

    Returns
    -------
    list of str
        The source flagged ancient, if any, then the frozen inputs it reads.
    """
    spec = freeze_spec(w.genome, w.name)
    source = [ancient(spec["source"])] if spec["source"] else []
    return [*source, *(frozen(w.genome, name) for name in spec["frozen"])]


# Input functions read the frozen lists to enumerate jobs, so a fresh inputs_dir
# needs this target before any other.
rule freeze:
    input:
        [
            frozen(genome, name)
            for genome in GENOMES
            for name in config["inputs"][genome]
        ],


rule freeze_input:
    input:
        input_sources,
    output:
        f"{INPUTS}/{{genome}}/{{name}}",
    params:
        kind=lambda w: freeze_spec(w.genome, w.name)["kind"],
        metadata=config["encode_metadata"],
        assembly=lambda w: config["genomes"][w.genome]["assembly"],
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
            --metadata {params.metadata:q} \
            --assembly {params.assembly:q} \
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
