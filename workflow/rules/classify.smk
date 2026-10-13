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
        """
