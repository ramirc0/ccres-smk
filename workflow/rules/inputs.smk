# Frozen inputs and ENCODE downloads.

# Derived inputs filter their sources. Every other input is a plain copy.
DERIVED = {
    "cCREs-Previous.bed": "previous-ccres",
    "TF-List.Cell-Type.txt": "cell-type-tf",
}


def input_sources(w):
    """Return the source files of a frozen input.

    Sources live in another user's tree and keep changing. `ancient` stops
    their mtimes from triggering reruns. The checksum guards their content.

    Parameters
    ----------
    w : snakemake.io.Wildcards
        Wildcards `genome` and `name`.

    Returns
    -------
    list of str
        Source paths, flagged ancient.
    """
    src = config["sources"][w.genome][w.name]
    return [ancient(s) for s in ([src] if isinstance(src, str) else src)]


# Input functions read the frozen lists to enumerate jobs, so a fresh inputs_dir
# needs this target before any other.
rule freeze:
    input:
        [frozen(genome, name) for genome in GENOMES for name in config["sources"][genome]],


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


# Only for files missing from the local mirror.
rule download_encode:
    output:
        f"{DOWNLOAD_DIR}/{{exp}}/{{acc}}.{{ext}}",
    log:
        f"{LOGDIR}/download_encode/{{exp}}/{{acc}}.{{ext}}.txt",
    benchmark:
        f"{BENCHDIR}/download_encode/{{exp}}/{{acc}}.{{ext}}.tsv",
    params:
        metadata=config["encode_metadata"],
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
