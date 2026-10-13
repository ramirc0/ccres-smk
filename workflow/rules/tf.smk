# Step 3: TF experiment list, peaks, rPeaks and TF clusters.


rule tf_list:
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
