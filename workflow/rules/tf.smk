# Step 3: TF experiment list, peaks, rPeaks and TF clusters.


# A checkpoint, because tf_peaks reads the filtered list to find its inputs.
checkpoint tf_list:
    input:
        lambda w: frozen(w.genome, "TF-List.All.txt"),
    output:
        f"{OUTDIR}/{{genome}}/TF/{{genome}}-TF-List.Filtered.txt",
    log:
        f"{LOGDIR}/tf_list/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/tf_list/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/tf_list.py \
            --input {input:q} \
            --output {output:q}
        """


def tf_peak_files(genome):
    """List the peak files of a genome's filtered TF experiments.

    Parameters
    ----------
    genome : str
        Genome key.

    Returns
    -------
    list of str
        narrowPeak path per experiment, excluded experiments left out.
    """
    excluded = config["genomes"][genome]["tf_excluded"]
    rows = read_rows(checkpoints.tf_list.get(genome=genome).output[0])
    return [encode_file(r[0], r[1], "bed.gz") for r in rows if r[0] not in excluded]


rule tf_peaks:
    input:
        tf_list=f"{OUTDIR}/{{genome}}/TF/{{genome}}-TF-List.Filtered.txt",
        peaks=lambda w: [ancient(p) for p in tf_peak_files(w.genome)],
    output:
        bed=f"{OUTDIR}/{{genome}}/TF/tmp.bed",
        summits=f"{OUTDIR}/{{genome}}/TF/tmp.summits",
    params:
        encode_root=ENCODE_ROOT,
        download_dir=DOWNLOAD_DIR,
        excluded=lambda w: config["genomes"][w.genome]["tf_excluded"],
    log:
        f"{LOGDIR}/tf_peaks/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/tf_peaks/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/tf_peaks.py \
            --tf-list {input.tf_list:q} \
            --encode-root {params.encode_root:q} \
            --download-dir {params.download_dir:q} \
            --exclude {params.excluded:q} \
            --bed {output.bed:q} \
            --summits {output.summits:q}
        """


rule tf_rpeaks:
    input:
        f"{OUTDIR}/{{genome}}/TF/tmp.bed",
    output:
        f"{OUTDIR}/{{genome}}/TF/tmp.rPeaks",
    log:
        f"{LOGDIR}/tf_rpeaks/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/tf_rpeaks/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/tf_rpeaks.py \
            --bed {input:q} \
            --output {output:q}
        """


rule tf_clusters:
    input:
        rpeaks=f"{OUTDIR}/{{genome}}/TF/tmp.rPeaks",
        summits=f"{OUTDIR}/{{genome}}/TF/tmp.summits",
        rdhs=f"{OUTDIR}/{{genome}}/rDHS/{{genome}}-rDHS-Filtered.bed",
    output:
        clusters=f"{OUTDIR}/{{genome}}/TF/{{genome}}-TF-Clusters.txt",
        no_overlap=f"{OUTDIR}/{{genome}}/TF/tmp.no-overlap",
        summary=f"{OUTDIR}/{{genome}}/TF/{{genome}}-rTFBS-Summary.txt",
        rtfbs=f"{OUTDIR}/{{genome}}/TF/{{genome}}-rTFBS.bed",
    log:
        f"{LOGDIR}/tf_clusters/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/tf_clusters/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/tf_clusters.py \
            --rpeaks {input.rpeaks:q} \
            --summits {input.summits:q} \
            --rdhs {input.rdhs:q} \
            --genome {wildcards.genome:q} \
            --clusters {output.clusters:q} \
            --no-overlap {output.no_overlap:q} \
            --summary {output.summary:q} \
            --rtfbs {output.rtfbs:q}
        """
