rule lofreq_preprocess:
    input:
        bam=rule.mark_duplicates.output.bam,
        ref=config["database"]["hg19"],
    output:
        bam=temp("align/lofreq/{sample}.lofreq_ready.bam")
    log:
        ".log/align/lofreq/{sample}.lofreq_prep.log"
    benchmark:
        ".log/align/lofreq/{sample}.lofreq_prep.bm"
    conda:
        config["conda"]["lofreq"]
    shell:
        "lofreq indelqual --dindel -f {input.ref} -o {output.bam} {input.bam} 2> {log}"


use rule samtools_sort_and_index as lofreq_bam_sort_and_index with:
    input:
        bam=rule.lofreq_preprocess.output.bam,
    output:
        bam="align/lofreq/{sample}.lofreq_ready_sorted.bam",
        bai="align/lofreq/{sample}.lofreq_ready_sorted.bam.bai",
    log:
        ".log/align/lofreq/{sample}.lofreq_sort_and_index.log"
    benchmark:
        ".log/align/lofreq/{sample}.lofreq_sort_and_index.bm"


rule lofreq_call:
    input:
        bam=rules.lofreq_bam_sort_and_index.output.bam,
        ref=config["database"]["hg19"]
    output:
        vcf="variants/lofreq/{sample}.vcf"
    log:
        ".log/variants/lofreq/{sample}.log"
    benchmark:
        ".log/variants/lofreq/{sample}.lofreq_call.bm"
    params:
        # min-bq 严格过滤低质量碱基; sig 显著性阈值; min-alt-bq 严格过滤低质量碱基
        extra="--min-bq 30 --sig 0.01 --min-alt-bq 30 --call-indels"
    conda:
        config["conda"]["lofreq"]
    threads: config["threads"]["low"]  # LoFreq 支持多线程通过 --parallel
    shell:
        "lofreq call-parallel --pp-threads {threads} -f {input.ref} -o {output.vcf} {params.extra} {input.bam} 2> {log}"


rule lofreq_filter:
    input:
        vcf=rule.lofreq_call.output.vcf,
    output:
        vcf="variants/lofreq/{sample}.filtered.vcf"
    log:
        ".log/variants/lofreq/{sample}.lofreq_filter.log"
    benchmark:
        ".log/variants/lofreq/{sample}.lofreq_filter.bm"
    conda:
        config["conda"]["lofreq"]
    params:
        # v 最小覆盖深度; a 最小等位基因频率
        extra="-v 5 -a 0.001"
    shell:
        "lofreq filter {input.vcf} -o {output.vcf} {params.extra} --call-indels 2> {log}"


rule variant_bgzip_and_index:
    input:
        vcf=rule.lofreq_filter.output.vcf,
    output:
        vcf="variants/lofreq/{sample}.filtered.vcf.gz",
        csi="variants/lofreq/{sample}.filtered.vcf.gz.csi",
    log:
        ".log/variants/lofreq/{sample}.variant_bgzip_and_index.log"
    benchmark:
        ".log/variants/lofreq/{sample}.variant_bgzip_and_index.bm"
    conda:
        config["conda"]["bcftools"]
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf} 2> {log}
        bcftools index {output.vcf} 2> {log}
        """
