rule bwa_mem:
    input:
        fq=rules.fastp.output.trimmed,
        ref=config["database"]["hg19"],
    output:
        "align/bwa/{sample}.sam",
    benchmark:
        ".log/align/bwa/{sample}.bwa_mem.bm"
    log:
        ".log/align/bwa/{sample}.bwa_mem.log",
    conda:
        config["conda"]["bwa"]
    threads: config["threads"]["medium"]
    params:
        extra=r"-M -Y -R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
    shell:
        "bwa mem -t {threads} {params.extra} {input.ref} {input.fq} -o {output} 2> {log}"


rule samtools_sort_and_index:
    input:
        rules.bwa_mem.output,
    output:
        bam="align/bwa/{sample}.bam",
        bai="align/bwa/{sample}.bam.bai",
    benchmark:
        ".log/align/bwa/{sample}.samtools_sort_and_index.bm"
    log:
        ".log/align/bwa/{sample}.samtools_sort_and_index.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["medium"]
    shell:
        """
        samtools sort -o {output.bam} {input} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


# BQSR 校正
rule recalibrate_base_qualities:
    input:
        bam=rules.samtools_sort_and_index.output.bam,
        bai=rules.samtools_sort_and_index.output.bai,
        ref=config["database"]["hg19"],
        gatk_dict=config["database"]["gatk_dict"],
        bed=f"{workflow.basedir}/assets/probeCov.predict.bed",
        known=[
            config["database"]["known_site_1000g"],
            config["database"]["known_site_dbsnp"],
            config["database"]["known_site_mills"],
        ],
        known_idx=[
            config["database"]["known_site_1000g_idx"],
            config["database"]["known_site_dbsnp_idx"],
            config["database"]["known_site_mills_idx"],
        ],
    output:
        recal_table="align/bqsr/{sample}-recal.grp",
    benchmark:
        ".log/align/bqsr/{sample}.recalibrate_base_qualities.bm"
    log:
        ".log/align/bqsr/{sample}.recalibrate_base_qualities.log",
    conda:
        config["conda"]["gatk"]
    shell:
        """
        gatk BaseRecalibrator \
            --input {input.bam} \
            --reference {input.ref} \
            --known-sites {input.known[0]} \
            --known-sites {input.known[1]} \
            --known-sites {input.known[2]} \
            --intervals {input.bed} \
            --output {output.recal_table} \
            2> {log}
        """


rule apply_base_quality_recalibration:
    input:
        bam=rules.samtools_sort_and_index.output.bam,
        bai=rules.samtools_sort_and_index.output.bai,
        ref=config["database"]["hg19"],
        gatk_dict=config["database"]["gatk_dict"],
        recal_table=rules.recalibrate_base_qualities.output.recal_table,
        bed=f"{workflow.basedir}/assets/probeCov.predict.bed",
    output:
        bam="align/bqsr/{sample}.recal.bam",
    benchmark:
        ".log/align/bqsr/{sample}.apply_base_quality_recalibration.bm"
    log:
        ".log/align/bqsr/{sample}.apply_base_quality_recalibration.log",
    conda:
        config["conda"]["gatk"]
    shell:
        """
        gatk ApplyBQSR \
            --input {input.bam} \
            --bqsr-recal-file {input.recal_table} \
            --reference {input.ref} \
            --intervals {input.bed} \
            --output {output.bam} \
            2> {log}
        """


use rule samtools_sort_and_index as recal_sort_and_index with:
    input:
        rules.apply_base_quality_recalibration.output.bam,
    output:
        bam="align/bqsr/{sample}.recal.sorted.bam",
        bai="align/bqsr/{sample}.recal.sorted.bam.bai",
    benchmark:
        ".log/align/bqsr/{sample}.recal_sort_and_index.bm"
    log:
        ".log/align/bqsr/{sample}.recal_sort_and_index.log",
