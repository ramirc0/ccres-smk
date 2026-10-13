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
