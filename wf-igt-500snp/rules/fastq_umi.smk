rule fastq_umi_pigz_decompress:
    input:
        rules.create_symlinks.output.fq1,
    output:
        temp("qc/umi/{sample}.1.fastq"),
    benchmark:
        ".log/qc/umi/{sample}.pigz_decompress.bm"
    log:
        ".log/qc/umi/{sample}.pigz_decompress.log",
    threads: config["threads"]["low"]
    shell:
        "pigz -d -c {input} > {output} 2> {log}"


use rule fastq_umi_pigz_decompress as fastq_umi_pigz_decompress_r2 with:
    input:
        rules.create_symlinks.output.fq2,
    output:
        temp("qc/umi/{sample}.2.fastq"),
    benchmark:
        ".log/qc/umi/{sample}.pigz_decompress.bm"
    log:
        ".log/qc/umi/{sample}.pigz_decompress.log",


"""

该 AI 生成脚本已经过人工核对, 确认无误
"""
rule fastq_umi_filter:
    input:
        r1=rules.fastq_umi_pigz_decompress.output,
        r2=rules.fastq_umi_pigz_decompress_r2.output,
        umi_file=f"{workflow.basedir}/assets/umi.txt",
    output:
        r1="qc/umi/{sample}.umi.1.fastq",
        r2="qc/umi/{sample}.umi.2.fastq",
        stats="qc/umi/{sample}.umi.stats.txt",
    message:
        "Filtering UMI for {input.r1} and {input.r2} to {output.r1} and {output.r2} with {input.umi_file} and {output.stats}"
    benchmark:
        ".log/qc/umi/{sample}.fastq_umi_filter.bm"
    log:
        ".log/qc/umi/{sample}.fastq_umi_filter.log",
    params:
        rust_umi_filter=f"{workflow.basedir}/tools/rust_umi_filter",
    shell:
        "{params.rust_umi_filter} {input.r1} {input.r2} {output.r1} {output.r2} {input.umi_file} {output.stats} 2> {log}"


rule fastq_umi_stats_summary:
    input:
        expand("qc/umi/{sample}.umi.stats.txt", sample=samples),
    output:
        "qc/umi/fastq_umi_stats_summary.tsv",
    benchmark:
        ".log/qc/umi/fastq_umi_stats_summary.bm"
    log:
        ".log/qc/umi/fastq_umi_stats_summary.log",
    script:
        "../scripts/fastq_umi_stats_summary.py"
