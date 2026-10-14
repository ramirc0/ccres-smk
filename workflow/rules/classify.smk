# Steps 6-8: cell type-agnostic cCREs, cell type-specific cCREs, PLS genes.


rule classify_ccres:
    input:
        anchors=f"{OUTDIR}/{{genome}}/{{genome}}-Anchors.bed",
        dnase=f"{OUTDIR}/{{genome}}/{{genome}}-DNase-maxZ.txt",
        h3k4me3=f"{OUTDIR}/{{genome}}/{{genome}}-H3K4me3-maxZ.txt",
        h3k27ac=f"{OUTDIR}/{{genome}}/{{genome}}-H3K27ac-maxZ.txt",
        ctcf=f"{OUTDIR}/{{genome}}/{{genome}}-CTCF-maxZ.txt",
        tss=lambda w: frozen(w.genome, "TSS.Basic.bed"),
        tss_4k=lambda w: frozen(w.genome, "TSS.Basic.4K.bed"),
        tf_clusters=f"{OUTDIR}/{{genome}}/TF/{{genome}}-TF-Clusters.txt",
        previous=lambda w: frozen(w.genome, "cCREs-Previous.bed"),
    output:
        f"{OUTDIR}/{{genome}}/{{genome}}-cCREs-Unfiltered.bed",
    log:
        f"{LOGDIR}/classify_ccres/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/classify_ccres/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/classify_ccres.py \
            --anchors {input.anchors:q} \
            --dnase {input.dnase:q} \
            --h3k4me3 {input.h3k4me3:q} \
            --h3k27ac {input.h3k27ac:q} \
            --ctcf {input.ctcf:q} \
            --tss {input.tss:q} \
            --tss-4k {input.tss_4k:q} \
            --tf-clusters {input.tf_clusters:q} \
            --previous {input.previous:q} \
            --genome {wildcards.genome:q} \
            --output {output:q}

        sync {output:q}
        """


import master_cell_list as mcl


@cache
def _master_rows(lists, tf_list, keys):
    """Return the master cell list rows of frozen lists, once per path set.

    Parameters
    ----------
    lists : tuple of str
        DNase, H3K4me3, H3K27ac and CTCF lists.
    tf_list : str
        Cell type TF experiment list.
    keys : str
        Biosample-Keys.txt.

    Returns
    -------
    list of list of str
        Master cell list rows, each with its biosample key last.
    """
    return mcl.master_rows(list(lists), tf_list, keys)


def master_rows(genome):
    """Return a genome's master cell list rows, built from its frozen lists.

    Parameters
    ----------
    genome : str
        Genome key.

    Returns
    -------
    list of list of str
        Master cell list rows, each with its biosample key last.
    """
    lists = tuple(frozen_list(genome, f"{mark}-List.txt") for mark in mcl.MARKS)
    return _master_rows(
        lists,
        frozen_list(genome, "TF-List.Cell-Type.txt"),
        frozen_list(genome, "Biosample-Keys.txt"),
    )


@cache
def tf_experiments(genome):
    """Group a genome's cell type TF experiments by biosample key.

    Parameters
    ----------
    genome : str
        Genome key.

    Returns
    -------
    dict of str to list of tuple of (str, str)
        (experiment, peak file) per key, tagging modifications aside.
    """
    keys = mcl.read_keys(frozen_list(genome, "Biosample-Keys.txt"))
    return mcl.tf_experiments(frozen_list(genome, "TF-List.Cell-Type.txt"), keys)


def cell_type_files(genome):
    """List the step 7 outputs of a genome, one per master cell list row.

    Parameters
    ----------
    genome : str
        Genome key.

    Returns
    -------
    list of str
        bgzipped BED and bigBed per biosample.
    """
    return [
        f"{OUTDIR}/{genome}/Cell-Type-Specific/Individual-Files/{mcl.output_name(r)}.{ext}"
        for r in master_rows(genome)
        for ext in ("bed.gz", "bigBed")
    ]


def cell_type_inputs(w):
    """Name the inputs of one biosample's classification.

    Parameters
    ----------
    w : snakemake.io.Wildcards
        Wildcards genome and cts.

    Returns
    -------
    dict of str to str or list of str
        Signal tables per present mark and, for the TF branch, the peak files
        of the biosample's TF experiments.
    """
    (row,) = [r for r in master_rows(w.genome) if mcl.output_name(r) == w.cts]
    inputs = {
        mark.lower(): f"{OUTDIR}/{w.genome}/signal-output/{mark}/{row[2 * i]}-{row[2 * i + 1]}.txt"
        for i, mark in enumerate(mcl.MARKS)
        if row[2 * i] != mcl.MISSING
    }
    if row[0] != mcl.MISSING and row[10] == "yes":
        inputs["tf_peaks"] = [
            ancient(encode_file(exp, peak, "bed.gz"))
            for exp, peak in tf_experiments(w.genome)[row[-1]]
        ]
    return inputs


def cell_type_flags(w, input):
    """Return the cell_type_ccres.py flags of one biosample.

    Parameters
    ----------
    w : snakemake.io.Wildcards
        Wildcards genome and cts.
    input : snakemake.io.InputFiles
        The job's named inputs.

    Returns
    -------
    list of str
        Signal flags, then the TF peak files of the TF branch.
    """
    flags = []
    for mark in mcl.MARKS:
        if mark.lower() in input.keys():
            flags += [f"--{mark.lower()}", input[mark.lower()]]
    if "tf_peaks" in input.keys():
        flags += ["--tf-peaks", *input.tf_peaks]
    return flags


rule master_cell_list:
    input:
        **{
            mark.lower(): (lambda w, mark=mark: frozen(w.genome, f"{mark}-List.txt"))
            for mark in mcl.MARKS
        },
        tf_list=lambda w: frozen(w.genome, "TF-List.Cell-Type.txt"),
        biosample_keys=lambda w: frozen(w.genome, "Biosample-Keys.txt"),
    output:
        f"{OUTDIR}/{{genome}}/Cell-Type-Specific/Master-Cell-List.txt",
    log:
        f"{LOGDIR}/master_cell_list/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/master_cell_list/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/master_cell_list.py \
            --dnase {input.dnase:q} \
            --h3k4me3 {input.h3k4me3:q} \
            --h3k27ac {input.h3k27ac:q} \
            --ctcf {input.ctcf:q} \
            --tf-list {input.tf_list:q} \
            --keys {input.biosample_keys:q} \
            --output {output:q}

        sync {output:q}
        """


rule cell_type_ccres:
    input:
        unpack(cell_type_inputs),
        master=f"{OUTDIR}/{{genome}}/{{genome}}-cCREs-Unfiltered.bed",
        tss=lambda w: frozen(w.genome, "TSS.Basic.bed"),
        tss_4k=lambda w: frozen(w.genome, "TSS.Basic.4K.bed"),
        autosql=lambda w: frozen(w.genome, "cres.as"),
        chrom_sizes=lambda w: frozen(w.genome, "chromInfo.txt"),
    output:
        bed=f"{OUTDIR}/{{genome}}/Cell-Type-Specific/Individual-Files/{{cts}}.bed.gz",
        tbi=f"{OUTDIR}/{{genome}}/Cell-Type-Specific/Individual-Files/{{cts}}.bed.gz.tbi",
        bigbed=f"{OUTDIR}/{{genome}}/Cell-Type-Specific/Individual-Files/{{cts}}.bigBed",
    params:
        flags=cell_type_flags,
    wildcard_constraints:
        cts=r"[A-Za-z0-9_]+",
    log:
        f"{LOGDIR}/cell_type_ccres/{{genome}}/{{cts}}.txt",
    benchmark:
        f"{BENCHDIR}/cell_type_ccres/{{genome}}/{{cts}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/cell_type_ccres.py \
            --master {input.master:q} \
            --tss {input.tss:q} \
            --tss-4k {input.tss_4k:q} \
            {params.flags:q} \
            --autosql {input.autosql:q} \
            --chrom-sizes {input.chrom_sizes:q} \
            --output {output.bed:q} \
            --bigbed {output.bigbed:q}

        sync {output:q}
        """


rule pls_genes:
    input:
        ccres=f"{OUTDIR}/{{genome}}/{{genome}}-cCREs-Unfiltered.bed",
        tss=lambda w: frozen(w.genome, "TSS.Basic.bed"),
    output:
        f"{OUTDIR}/{{genome}}/{{genome}}-PLS-Gene-List.txt",
    log:
        f"{LOGDIR}/pls_genes/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/pls_genes/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/pls_genes.py \
            --ccres {input.ccres:q} \
            --tss {input.tss:q} \
            --output {output:q}

        sync {output:q}
        """
