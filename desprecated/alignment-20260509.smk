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


# ! 重复标记为1024(0x400)，在 IGV 中自动被忽略. 所以会显示深度与实际不符. 看是否要删除这一步.
rule mark_duplicates:
    input:
        bam=rules.samtools_sort_and_index.output.bam,
    output:
        bam="align/markdup/{sample}.md.bam",  # 标记后的中间文件
        metrics="align/markdup/{sample}.md_metrics.txt",
    log:
        ".log/align/markdup/{sample}.mark_duplicates.log",
    benchmark:
        ".log/align/markdup/{sample}.mark_duplicates.bm"
    conda:
        config["conda"]["gatk"]
    threads: config["threads"]["low"]
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --REMOVE_DUPLICATES false \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT > {log} 2>&1
        """
