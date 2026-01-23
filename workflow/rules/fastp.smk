rule fastp:
    input:
        sample=[
            rules.create_symlinks.output.fq1,
            rules.create_symlinks.output.fq2,
        ],
    output:
        trimmed=["qc/fastp/{sample}.1.fastq.gz", "qc/fastp/{sample}.2.fastq.gz"],
        html="qc/fastp/{sample}.html",
        json="qc/fastp/{sample}.json",
    log:
        ".log/qc/fastp/{sample}.fastp.log",
    benchmark:
        ".log/qc/fastp/{sample}.fastp.bm"
    conda:
        config["conda"]["fastp"]
    threads: config["threads"]["low"]
    params:
        extra="-q 15 -u 40 -l 25 --cut_right --cut_window_size 20 --cut_mean_quality 30 --correction",
    shell:
        """
        fastp -w {threads} {params.extra} \
            -i {input.sample[0]} -I {input.sample[1]} \
            -o {output.trimmed[0]} -O {output.trimmed[1]} \
            -h {output.html} -j {output.json} 2> {log}
        """


rule fq_stats_summary:
    input:
        expand("fastp/{sample}.json", sample=samples),
    output:
        "qc/fastp/fq_summary.tsv",
    benchmark:
        ".log/qc/fastp/fq_all_samples_qc.bm"
    log:
        ".log/qc/fastp/fq_all_samples_qc.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/fq_all_samples_qc.py"
