# Anchors, step 4 (z-scores) and step 5 (max z-scores).

# bp added to each side of an anchor before averaging, per Retrieve-Signal.sh.
WIDTH = {"DNase": 0, "CTCF": 0, "ATAC": 0, "H3K4me3": 500, "H3K27ac": 500}


def mark_signal_files(genome, mark):
    """List the step 4 outputs of one mark.

    Parameters
    ----------
    genome : str
        Genome key.
    mark : str
        Mark, a key of `WIDTH`.

    Returns
    -------
    list of str
        One z-score table per row of the mark list. Rows whose file is NA
        have none.
    """
    rows = read_rows(frozen_list(genome, f"{mark}-List.txt"))
    return [
        f"{OUTDIR}/{genome}/signal-output/{mark}/{r[0]}-{r[1]}.parquet"
        for r in rows
        if r[1] != "NA"
    ]


rule anchors:
    input:
        rdhs=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-rDHS-Filtered.bed",
        rtfbs=f"{OUTDIR}/{{genome}}/TF/{{genome}}-rTFBS.bed",
    output:
        f"{OUTDIR}/{{genome}}/{{genome}}-Anchors.parquet",
    log:
        f"{LOGDIR}/anchors/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/anchors/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/anchors.py \
            --rdhs {input.rdhs:q} \
            --rtfbs {input.rtfbs:q} \
            --output {output:q}

        sync {output:q}
        """


rule extend_anchors:
    input:
        f"{OUTDIR}/{{genome}}/{{genome}}-Anchors.parquet",
    output:
        f"{OUTDIR}/{{genome}}/signal-output/Anchors-{{width}}bp.parquet",
    wildcard_constraints:
        width=r"\d+",
    log:
        f"{LOGDIR}/extend_anchors/{{genome}}/{{width}}.txt",
    benchmark:
        f"{BENCHDIR}/extend_anchors/{{genome}}/{{width}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/extend_anchors.py \
            --anchors {input:q} \
            --width {wildcards.width:q} \
            --output {output:q}

        sync {output:q}
        """


rule signal_zscore:
    input:
        anchors=lambda w: f"{OUTDIR}/{w.genome}/signal-output/Anchors-{WIDTH[w.mark]}bp.parquet",
        bigwig=lambda w: ancient(encode_file(w.exp, w.signal, "bigWig")),
    output:
        f"{OUTDIR}/{{genome}}/signal-output/{{mark}}/{{exp}}-{{signal}}.parquet",
    wildcard_constraints:
        signal=r"ENCFF\d{3}[A-Z]{3}",
    log:
        f"{LOGDIR}/signal_zscore/{{genome}}/{{mark}}/{{exp}}-{{signal}}.txt",
    benchmark:
        f"{BENCHDIR}/signal_zscore/{{genome}}/{{mark}}/{{exp}}-{{signal}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/signal_zscore.py \
            --anchors {input.anchors:q} \
            --bigwig {input.bigwig:q} \
            --output {output:q}

        sync {output:q}
        """


rule max_zscore:
    input:
        mark_list=lambda w: frozen(w.genome, f"{w.mark}-List.txt"),
        signal=lambda w: mark_signal_files(w.genome, w.mark),
    output:
        f"{OUTDIR}/{{genome}}/{{genome}}-{{mark}}-maxZ.parquet",
    params:
        signal=lambda w: f"{OUTDIR}/{w.genome}/signal-output/{w.mark}",
    threads: 8
    log:
        f"{LOGDIR}/max_zscore/{{genome}}/{{mark}}.txt",
    benchmark:
        f"{BENCHDIR}/max_zscore/{{genome}}/{{mark}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/max_zscore.py \
            --mark-list {input.mark_list:q} \
            --signal {params.signal:q} \
            --threads {threads} \
            --output {output:q}

        sync {output:q}
        """
