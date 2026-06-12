rule bam_clipOverlap:
    input:
        bam=rules.recal_sort_and_index.output.bam,
        bai=rules.recal_sort_and_index.output.bai,
    output:
        temp("allele/{sample}.clipOverlap.bam"),
    benchmark:
        ".log/allele/{sample}.bam_clipOverlap.bm"
    log:
        ".log/allele/{sample}.bam_clipOverlap.log",
    conda:
        config["conda"]["bamutils"]
    shell:
        "bam clipOverlap --in {input.bam} --out {output} --stats 2> {log}"


use rule samtools_sort_and_index as bam_clipOverlap_sort_and_index with:
    input:
        rules.bam_clipOverlap.output,
    output:
        bam="allele/{sample}.clipOverlap.sorted.bam",
        bai="allele/{sample}.clipOverlap.sorted.bam.bai",
    benchmark:
        ".log/allele/{sample}.clipOverlap.sorted.bam.bm"
    log:
        ".log/allele/{sample}.clipOverlap.sorted.bam.log",


rule allele_mpileup:
    input:
        bam=rules.bam_clipOverlap_sort_and_index.output.bam,
        bai=rules.bam_clipOverlap_sort_and_index.output.bai,
        reference=config["database"]["hg19"],
        positions=f"{workflow.basedir}/assets/loci.position.txt",
    output:
        mpileup="allele/{sample}.clipOverlap.mpileup",
        tmpdir=temp(directory("allele/tmp-{sample}")),
    benchmark:
        ".log/allele/{sample}.clipOverlap.mpileup.bm"
    log:
        ".log/allele/{sample}.clipOverlap.mpileup.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["medium"]
    params:
        "--max-depth 500000 --min-BQ 30",
    shell:
        """
        mkdir -p {output.tmpdir}
        cat {input.positions} | while read position; do
            samtools mpileup \
            --fasta-ref {input.reference} \
            {params} \
            -r $position {input.bam} > {output.tmpdir}/$position.mpileup
        done | parallel -j {threads} 2> {log}
        cat {output.tmpdir}/*.mpileup > {output.mpileup}
        """


rule parse_bases:
    input:
        mpileup=rules.allele_mpileup.output.mpileup,
    output:
        base_stat="allele/{sample}.clipOverlap.base_stat.tsv",
    benchmark:
        ".log/allele/{sample}.parse_bases.bm"
    log:
        ".log/allele/{sample}.parse_bases.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/parse_bases.py"


rule calc_af:
    input:
        vcf=f"{workflow.basedir}/assets/loci.vcf",
        base_stat=rules.parse_bases.output.base_stat,
    output:
        alt_freq="allele/{sample}.alt_freq.tsv",
    benchmark:
        ".log/allele/{sample}.calc_af.bm"
    log:
        ".log/allele/{sample}.calc_af.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/calc_af.py"


rule calc_ddcfDNA_pct:
    input:
        expand("allele/{sample}.alt_freq.tsv", sample=samples),
    output:
        tsv="upload/all.dd_cfDNA_pct.tsv",
        xlsx= "upload/all.dd_cfDNA_pct.xlsx",
    benchmark:
        ".log/upload/all.dd_cfDNA_pct.bm"
    log:
        ".log/upload/all.dd_cfDNA_pct.log",
    conda:
        config["conda"]["python"]
    params:
        method="both",  # ['caredx', 'gmm', 'both']
        multiplier=2.11,  # 亲缘关系校正乘数 M (default: 2.11, 无关供受者)
        min_depth=1000,  # 最小深度 (default: 1000)
        het_threshold=0.35,  # HET 阈值 (default: 0.35)
        outlier_fraction=0.05,  # 异常值分数 (default: 0.05)
        max_contam = 0.01  # 最大 CONTAM_FREQ 过滤 (default: 0.01)
    script:
        "../scripts/calc_ddcfdfna.py"
