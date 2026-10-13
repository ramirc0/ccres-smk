# Steps 0-1: DHSs of each DNase experiment, then their mean signal.


def processed_dhs(genome):
    """List the step 1 outputs of a genome, one per hotspot list row.

    Parameters
    ----------
    genome : str
        Genome key.

    Returns
    -------
    list of str
        Processed DHS BED per row.
    """
    rows = read_rows(frozen(genome, "Hotspot-List.txt"))
    return [
        f"{OUTDIR}/{genome}/rDHS/Processed-DHSs/{r[2]}.{i}.bed"
        for i, r in enumerate(rows, start=1)
    ]


rule call_dhs:
    input:
        lambda w: ancient(
            encode_file(hotspot_experiment(w.genome, w.acc), w.acc, "bed.gz")
        ),
    output:
        dhs=f"{OUTDIR}/{{genome}}/rDHS/DHSs/{{acc}}.DHSs.bed",
        excluded=f"{OUTDIR}/{{genome}}/rDHS/DHSs/{{acc}}.Excluded.bed",
    params:
        min_p=lambda w: config["genomes"][w.genome]["min_p"],
    log:
        f"{LOGDIR}/call_dhs/{{genome}}/{{acc}}.txt",
    benchmark:
        f"{BENCHDIR}/call_dhs/{{genome}}/{{acc}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/call_dhs.py \
            --hotspots {input:q} \
            --min-p {params.min_p:q} \
            --dhs {output.dhs:q} \
            --excluded {output.excluded:q}
        """


# {acc} is the signal file of hotspot list row {row}. The original named its
# outputs output.<acc>.<row>, and step 2 concatenates them in that name order.
rule process_dhs:
    input:
        dhs=lambda w: f"{OUTDIR}/{w.genome}/rDHS/DHSs/{hotspot_row(w.genome, w.row)[1]}.DHSs.bed",
        bigwig=lambda w: ancient(
            encode_file(hotspot_row(w.genome, w.row)[0], w.acc, "bigWig")
        ),
    output:
        f"{OUTDIR}/{{genome}}/rDHS/Processed-DHSs/{{acc}}.{{row}}.bed",
    params:
        peak=lambda w: hotspot_row(w.genome, w.row)[1],
    log:
        f"{LOGDIR}/process_dhs/{{genome}}/{{acc}}.{{row}}.txt",
    benchmark:
        f"{BENCHDIR}/process_dhs/{{genome}}/{{acc}}.{{row}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/process_dhs.py \
            --dhs {input.dhs:q} \
            --peak {params.peak:q} \
            --bigwig {input.bigwig:q} \
            --output {output:q}
        """


rule combine_dhs:
    input:
        hotspots=lambda w: frozen(w.genome, "Hotspot-List.txt"),
        processed=lambda w: processed_dhs(w.genome),
    output:
        all=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-DHS-All.bed",
        filtered=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-DHS-Filtered.bed",
    params:
        processed=lambda w: f"{OUTDIR}/{w.genome}/rDHS/Processed-DHSs",
        cutoff=lambda w: config["genomes"][w.genome]["dhs_signal_cutoff"],
    log:
        f"{LOGDIR}/combine_dhs/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/combine_dhs/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/combine_dhs.py \
            --hotspots {input.hotspots:q} \
            --processed {params.processed:q} \
            --cutoff {params.cutoff:q} \
            --all {output.all:q} \
            --filtered {output.filtered:q}
        """


rule create_rdhs:
    input:
        filtered=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-DHS-Filtered.bed",
        previous=lambda w: frozen(w.genome, "rDHS-Previous.bed.gz"),
        cdhs=lambda w: frozen(w.genome, "cDHSs.bed"),
    output:
        prestam=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-rDHS-PreStam.txt",
        complete=f"{OUTDIR}/{{genome}}/rDHS/completeCoverage.bed",
        none=f"{OUTDIR}/{{genome}}/rDHS/NoCoverage.bed",
        summary=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-rDHS-Filtered-Summary.txt",
        rdhs=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-rDHS-Filtered.bed",
    log:
        f"{LOGDIR}/create_rdhs/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/create_rdhs/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/create_rdhs.py \
            --filtered {input.filtered:q} \
            --previous {input.previous:q} \
            --cdhs {input.cdhs:q} \
            --genome {wildcards.genome:q} \
            --prestam {output.prestam:q} \
            --complete {output.complete:q} \
            --none {output.none:q} \
            --summary {output.summary:q} \
            --rdhs {output.rdhs:q}
        """


# Reports the data-derived values of the hardcoded constants. The run keeps the
# release values.
rule derive_constants:
    input:
        hotspots=lambda w: frozen(w.genome, "Hotspot-List.txt"),
        peaks=lambda w: [
            ancient(encode_file(r[0], r[1], "bed.gz"))
            for r in read_rows(frozen(w.genome, "Hotspot-List.txt"))
        ],
        dhs_all=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-DHS-All.bed",
    output:
        f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-Constants.tsv",
    params:
        encode_root=ENCODE_ROOT,
        download_dir=DOWNLOAD_DIR,
        min_p=lambda w: config["genomes"][w.genome]["min_p"],
        cutoff=lambda w: config["genomes"][w.genome]["dhs_signal_cutoff"],
    threads: 8
    log:
        f"{LOGDIR}/derive_constants/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/derive_constants/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/derive_constants.py \
            --hotspots {input.hotspots:q} \
            --encode-root {params.encode_root:q} \
            --download-dir {params.download_dir:q} \
            --dhs-all {input.dhs_all:q} \
            --min-p {params.min_p:q} \
            --cutoff {params.cutoff:q} \
            --threads {threads:q} \
            --output {output:q}
        """
