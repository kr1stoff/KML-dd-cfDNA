rule samtools_stats:
    input:
        bam=rules.bwa_mem_raw.output.bam,
        bed=config["database"]["region"],
        idx=rules.bwa_mem_raw.output.bai,
    output:
        "align/stats/{sample}.bam.target.stat",
    benchmark:
        ".log/align/stats/{sample}.samtools_stats.bm"
    log:
        ".log/align/stats/{sample}.samtools_stats.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["low"]
    shell:
        "samtools stats --threads {threads} --target-regions {input.bed} {input.bam} > {output} 2> {log}"


rule samtools_stats_all:
    input:
        rules.bwa_mem_raw.output.bam,
    output:
        "align/stats/{sample}.bam.stat",
    benchmark:
        ".log/align/stats/{sample}.samtools_stats_all.bm"
    log:
        ".log/align/stats/{sample}.samtools_stats_all.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["low"]
    shell:
        "samtools stats --threads {threads} {input} > {output} 2> {log}"


rule samtools_depth:
    input:
        bam=rules.bwa_mem_raw.output.bam,
        bed=config["database"]["region"],
    output:
        "align/stats/{sample}.bam.target.depth",
    benchmark:
        ".log/align/stats/{sample}.samtools_depth.bm"
    log:
        ".log/align/stats/{sample}.samtools_depth.log",
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools depth -a -b {input.bed} {input.bam} > {output} 2> {log}"


rule bam_stats:
    input:
        rules.samtools_stats_all.output,
        rules.samtools_stats.output,
        rules.samtools_depth.output,
    output:
        "align/stats/{sample}.stats.csv",
    benchmark:
        ".log/align/stats/{sample}.bam_stats.bm"
    log:
        ".log/align/stats/{sample}.bam_stats.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats.py"


rule bam_stats_summary:
    input:
        expand("align/stats/{sample}.stats.csv", sample=samples),
    output:
        "align/stats/bam_summary.tsv",
    benchmark:
        ".log/align/stats/bam_stats_summary.bm"
    log:
        ".log/align/stats/bam_stats_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats_summary.py"
