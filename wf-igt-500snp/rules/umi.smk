rule picard_fastq_to_sam:
    input:
        fq1=rules.fastp.output.trimmed[0],
        fq2=rules.fastp.output.trimmed[1],
        tmp_dir=config["resources"]["tmp_dir"],
    output:
        temp("align/umi/{sample}.unmapped.bam"),
    benchmark:
        ".log/align/umi/{sample}.picard_fastq_to_sam.bm"
    log:
        ".log/align/umi/{sample}.picard_fastq_to_sam.log",
    conda:
        config["conda"]["picard"]
    params:
        java_mem=config["resources"]["java_mem"],
    shell:
        "picard FastqToSam {params.java_mem} F1={input.fq1} F2={input.fq2} O={output} SAMPLE_NAME={wildcards.sample} TMP_DIR={input.tmp_dir} 2> {log}"


rule fgbio_extract_umis:
    input:
        bam=rules.picard_fastq_to_sam.output,
    output:
        temp("align/umi/{sample}.unmapped.withUMI.bam"),
    benchmark:
        ".log/align/umi/{sample}.fgbio_extract_umis.bm"
    log:
        ".log/align/umi/{sample}.fgbio_extract_umis.log",
    conda:
        config["conda"]["fgbio"]
    params:
        java_mem=config["resources"]["java_mem"],
        extra="--read-structure=6M144T 6M144T --molecular-index-tags=ZA ZB --single-tag=RX",
    shell:
        "fgbio ExtractUmisFromBam {params.java_mem} --input={input.bam} --output={output} {params.extra} 2> {log}"


rule picard_sam_to_fastq:
    input:
        bam=rules.fgbio_extract_umis.output,
        tmp_dir=config["resources"]["tmp_dir"],
    output:
        temp("align/umi/{sample}.SamToFastq"),
    benchmark:
        ".log/align/umi/{sample}.picard_sam_to_fastq.bm"
    log:
        ".log/align/umi/{sample}.picard_sam_to_fastq.log",
    conda:
        config["conda"]["picard"]
    params:
        java_mem=config["resources"]["java_mem"],
        extra="INTERLEAVE=true",
    shell:
        "picard SamToFastq {params.java_mem} I={input.bam} F={output} {params.extra} TMP_DIR={input.tmp_dir} 2> {log}"


rule bwa_mem_raw:
    input:
        fq=rules.picard_sam_to_fastq.output,
        ref=config["database"]["hg19"],
    output:
        bam="align/umi/{sample}.mapped.bam",
        bai="align/umi/{sample}.mapped.bam.bai",
    benchmark:
        ".log/align/umi/{sample}.bwa_mem_raw.bm"
    log:
        ".log/align/umi/{sample}.bwa_mem_raw.log",
    conda:
        config["conda"]["bwa"]
    threads: config["threads"]["medium"]
    shell:
        """
        bwa mem -p -t {threads} {input.ref} {input.fq} 2> {log} | \
            samtools view -@ {threads} -hbS - 2>> {log} | \
            samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        """


rule picard_merge_bam_alignment:
    input:
        unmapped=rules.fgbio_extract_umis.output,
        aligned=rules.bwa_mem_raw.output.bam,
        ref=config["database"]["hg19"],
        tmp_dir=config["resources"]["tmp_dir"],
    output:
        temp("align/umi/{sample}.mapped.withUMI.bam"),
    benchmark:
        ".log/align/umi/{sample}.picard_merge_bam_alignment.bm"
    log:
        ".log/align/umi/{sample}.picard_merge_bam_alignment.log",
    conda:
        config["conda"]["picard"]
    params:
        java_mem=config["resources"]["java_mem"],
        extra="SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true",
    shell:
        "picard MergeBamAlignment {params.java_mem} UNMAPPED={input.unmapped} ALIGNED={input.aligned} O={output} R={input.ref} {params.extra} TMP_DIR={input.tmp_dir} 2> {log}"


rule fgbio_group_reads_by_umi:
    input:
        bam=rules.picard_merge_bam_alignment.output,
    output:
        temp("align/umi/{sample}.grouped.bam"),
    benchmark:
        ".log/align/umi/{sample}.fgbio_group_reads_by_umi.bm"
    log:
        ".log/align/umi/{sample}.fgbio_group_reads_by_umi.log",
    conda:
        config["conda"]["fgbio"]
    params:
        java_mem=config["resources"]["java_mem"],
        extra="--strategy=paired --edits=1 --min-map-q=20",
    threads: config["threads"]["medium"]
    shell:
        "fgbio GroupReadsByUmi {params.java_mem} --threads={threads} --input={input.bam} --output={output} {params.extra} 2> {log}"


rule fgbio_call_duplex_consensus_reads:
    input:
        bam=rules.fgbio_group_reads_by_umi.output,
        tmp_dir=config["resources"]["tmp_dir"],
    output:
        temp("align/umi/{sample}.consensus.unmapped.bam"),
    benchmark:
        ".log/align/umi/{sample}.fgbio_call_duplex_consensus_reads.bm"
    log:
        ".log/align/umi/{sample}.fgbio_call_duplex_consensus_reads.log",
    conda:
        config["conda"]["fgbio"]
    params:
        java_mem=config["resources"]["java_mem"],
        extra="--error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=30",
    threads: config["threads"]["medium"]
    shell:
        "fgbio CallDuplexConsensusReads {params.java_mem} --threads={threads} --input={input.bam} --output={output} {params.extra} --tmp-dir={input.tmp_dir} 2> {log}"


rule fgbio_filter_consensus_reads:
    input:
        bam=rules.fgbio_call_duplex_consensus_reads.output,
        ref=config["database"]["hg19"],
        tmp_dir=config["resources"]["tmp_dir"],
    output:
        temp("align/umi/{sample}.consensus.filtered.unmapped.bam"),
    benchmark:
        ".log/align/umi/{sample}.fgbio_filter_consensus_reads.bm"
    log:
        ".log/align/umi/{sample}.fgbio_filter_consensus_reads.log",
    conda:
        config["conda"]["fgbio"]
    params:
        java_mem=config["resources"]["java_mem"],
        extra="--min-reads=2 1 1 --max-read-error-rate=0.05 --max-base-error-rate=0.1 --min-base-quality=50 --max-no-calls=0.05",
    shell:
        "fgbio FilterConsensusReads {params.java_mem} --input={input.bam} --output={output} --ref={input.ref} {params.extra} --tmp-dir={input.tmp_dir} 2> {log}"


use rule picard_sam_to_fastq as picard_sam_to_fastq_filter_consensus_reads with:
    input:
        rules.fgbio_filter_consensus_reads.output,
    output:
        temp("align/umi/{sample}.consensus.filtered.unmapped.Fastq"),
    benchmark:
        ".log/align/umi/{sample}.picard_sam_to_fastq_filter_consensus_reads.bm"
    log:
        ".log/align/umi/{sample}.picard_sam_to_fastq_filter_consensus_reads.log",


use rule bwa_mem_raw as bwa_mem_filter_consensus_reads with:
    input:
        fq=rules.picard_sam_to_fastq_filter_consensus_reads.output,
        ref=config["database"]["hg19"],
    output:
        bam=temp("align/umi/{sample}.consensus.filtered.mapped.bam"),
        bai=temp("align/umi/{sample}.consensus.filtered.mapped.bam.bai"),
    benchmark:
        ".log/align/umi/{sample}.bwa_mem_filter_consensus_reads.bm"
    log:
        ".log/align/umi/{sample}.bwa_mem_filter_consensus_reads.log",


rule picard_sort_sam_filter_consensus_mapped:
    input:
        bam=rules.bwa_mem_filter_consensus_reads.output.bam,
        tmp_dir=config["resources"]["tmp_dir"],
    output:
        temp("align/umi/{sample}.consensus.filtered.mapped.sorted.bam"),
    benchmark:
        ".log/align/umi/{sample}.picard_sort_sam_filter_consensus_mapped.bm"
    log:
        ".log/align/umi/{sample}.picard_sort_sam_filter_consensus_mapped.log",
    conda:
        config["conda"]["picard"]
    params:
        java_mem=config["resources"]["java_mem"],
        extra="SORT_ORDER=queryname",
    shell:
        "picard SortSam {params.java_mem} INPUT={input.bam} OUTPUT={output} {params.extra} TMP_DIR={input.tmp_dir} 2> {log}"


use rule picard_sort_sam_filter_consensus_mapped as picard_sort_sam_filter_consensus_unmapped with:
    input:
        rules.fgbio_filter_consensus_reads.output,
    output:
        temp("align/umi/{sample}.consensus.filtered.unmapped.sorted.bam"),
    benchmark:
        ".log/align/umi/{sample}.picard_sort_sam_filter_consensus_unmapped.bm"
    log:
        ".log/align/umi/{sample}.picard_sort_sam_filter_consensus_unmapped.log",


use rule picard_merge_bam_alignment as picard_merge_bam_alignment_filter_consensus_reads with:
    input:
        unmapped=rules.picard_sort_sam_filter_consensus_unmapped.output,
        aligned=rules.picard_sort_sam_filter_consensus_mapped.output,
        ref=config["database"]["hg19"],
    output:
        temp("align/umi/{sample}.mapped.withUMI.consensus.filtered.bam"),
    benchmark:
        ".log/align/umi/{sample}.picard_merge_bam_alignment_filter_consensus_reads.bm"
    log:
        ".log/align/umi/{sample}.picard_merge_bam_alignment_filter_consensus_reads.log",


rule fgbio_clip_bam:
    input:
        bam=rules.picard_merge_bam_alignment_filter_consensus_reads.output,
        ref=config["database"]["hg19"],
        tmp_dir=config["resources"]["tmp_dir"],
    output:
        "align/umi/{sample}.consensus.filtered.mapped.withUMI.clipped.bam",
    benchmark:
        ".log/align/umi/{sample}.fgbio_clip_bam.bm"
    log:
        ".log/align/umi/{sample}.fgbio_clip_bam.log",
    conda:
        config["conda"]["fgbio"]
    params:
        samtools="-n -u",
        java_mem=config["resources"]["java_mem"],
        fgbio="--clipping-mode=Soft --clip-overlapping-reads=true",
    shell:
        """
        samtools sort {params.samtools} {input.bam} 2> {log} | \
            fgbio ClipBam {params.java_mem} {params.fgbio} --input=/dev/stdin --output={output} --ref={input.ref} --tmp-dir={input.tmp_dir} 2>> {log}
        """


rule sort_umi_bam:
    input:
        rules.fgbio_clip_bam.output,
    output:
        bam="align/umi/{sample}.consensus.filtered.mapped.withUMI.clipped.sorted.bam",
        bai="align/umi/{sample}.consensus.filtered.mapped.withUMI.clipped.sorted.bam.bai",
    benchmark:
        ".log/align/umi/{sample}.sort_umi_bam.bm"
    log:
        ".log/align/umi/{sample}.sort_umi_bam.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["medium"]
    shell:
        """
        samtools sort --threads {threads} {input} -o {output.bam}
        samtools index {output.bam}
        """
