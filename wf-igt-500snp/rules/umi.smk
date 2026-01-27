rule picard_fastq_to_sam:
    input:
        rules.fastp.output.trimmed[0],
        rules.fastp.output.trimmed[1],
    output:
        "align/umi/{sample}.unmapped.bam",
    benchmark:
        ".log/align/umi/{sample}.picard_fastq_to_sam.bm"
    log:
        ".log/align/umi/{sample}.picard_fastq_to_sam.log",
    conda:
        config["conda"]["picard"]
    params:
        java_mem=config["custom"]["java_mem"],
    shell:
        "picard FastqToSam {params.java_mem} F1={input[0]} F2={input[1]} O={output} SAMPLE_NAME={wildcards.sample}"


rule fgbio_extract_umis:
    input:
        rules.picard_fastq_to_sam.output,
    output:
        "align/umi/{sample}.unmapped.withUMI.bam",
    benchmark:
        ".log/align/umi/{sample}.fgbio_extract_umis.bm"
    log:
        ".log/align/umi/{sample}.fgbio_extract_umis.log",
    conda:
        config["conda"]["fgbio"]
    params:
        java_mem=config["custom"]["java_mem"],
        extra: "--read-structure=6M144T 6M144T --molecular-index-tags=ZA ZB --single-tag=RX"
    shell:
        "/home/mengxf/miniforge3/envs/basic2/bin/fgbio ExtractUmisFromBam {params.java_mem} --input={input} --output={output} {params.extra}"


rule picard_sam_to_fastq:
    input:
        rules.fgbio_extract_umis.output,
    output:
        "align/umi/{sample}.SamToFastq",
    benchmark:
        ".log/align/umi/{sample}.picard_sam_to_fastq.bm"
    log:
        ".log/align/umi/{sample}.picard_sam_to_fastq.log",
    conda:
        config["conda"]["picard"]
    params:
        java_mem=config["custom"]["java_mem"],
        extra: "INTERLEAVE=true"
    shell:
        "picard SamToFastq {params.java_mem} I={input} F={output} {params.extra}"


rule bwa_mem_raw:
    input:
        rules.picard_sam_to_fastq.output,
    output:
        "align/umi/{sample}.mapped.bam",
    benchmark:
        ".log/align/umi/{sample}.bwa_mem.bm"
    log:
        ".log/align/umi/{sample}.bwa_mem.log",
    conda:
        config["conda"]["bwa"]
    threads:
        config["threads"]["medium"]
    params:

    shell:
        "bwa mem -p -t {threads} /data/mengxf/Database/reference/hg19/hg19.fa {input} > {output}"

# bwa mem -p -t 32 /data/mengxf/Database/reference/hg19/hg19.fa ZQJ12.SamToFastq > ZQJ12.mapped.bam

# /home/mengxf/miniforge3/envs/basic2/bin/picard MergeBamAlignment \
#     UNMAPPED=ZQJ12.unmapped.withUMI.bam \
#     ALIGNED=ZQJ12.mapped.bam \
#     O=ZQJ12.mapped.withUMI.bam \
#     R=/data/mengxf/Database/reference/hg19/hg19.fa \
#     SO=coordinate \
#     ALIGNER_PROPER_PAIR_FLAGS=true \
#     MAX_GAPS=-1 \
#     ORIENTATIONS=FR \
#     VALIDATION_STRINGENCY=SILENT \
#     CREATE_INDEX=true