# QC
rule fastqc:
    input:
        rules.create_symlinks.output.fq1,
        rules.create_symlinks.output.fq2,
    output:
        directory("qc/fastqc/{sample}"),
    benchmark:
        ".log/qc/fastqc/{sample}.bm"
    log:
        ".log/qc/fastqc/{sample}.log",
    conda:
        config["conda"]["fastqc"]
    threads: config["threads"]["low"]
    shell:
        "mkdir {output} && fastqc {input} -o {output} -t {threads} --extract &> {log}"


rule multiqc:
    input:
        expand("qc/fastqc/{sample}", sample=samples),
    output:
        directory("qc/multiqc"),
    benchmark:
        ".log/qc/multiqc/multiqc.bm"
    log:
        ".log/qc/multiqc/multiqc.log",
    conda:
        config["conda"]["multiqc"]
    shell:
        "multiqc {input} --outdir {output} 2> {log}"


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
        # UMI 数据尽量少剪切, 否则会导致 UMI 丢失
        extra="",
    shell:
        """
        fastp -w {threads} {params.extra} \
            -i {input.sample[0]} -I {input.sample[1]} \
            -o {output.trimmed[0]} -O {output.trimmed[1]} \
            -h {output.html} -j {output.json} 2> {log}
        """


rule fq_stats_summary:
    input:
        expand("qc/fastp/{sample}.json", sample=samples),
    output:
        "qc/fastp/fq_summary.tsv",
    benchmark:
        ".log/qc/fastp/fq_all_samples_qc.bm"
    log:
        ".log/qc/fastp/fq_all_samples_qc.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/fq_stats_summary.py"
