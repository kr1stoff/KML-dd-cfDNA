rule bwa_mem_hg19_stat:
    input:
        rules.bwa_mem_hg19.output,
    output:
        "align/stats/{sample}.bam.stat",
    benchmark:
        ".log/align/stats/{sample}.bwa_mem_hg19_stat.bm"
    log:
        ".log/align/stats/{sample}.bwa_mem_hg19_stat.log",
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools stat {output.bam} | grep ^SN | cut -f 2- > {output} 2>> {log}"


rule filter_bam:
    input:
        rules.bwa_mem_hg19.output,
    output:
        bam=temp("align/filter/{sample}.filter.bam"),
        stat="align/stats/{sample}.filter.bam.stat",
    benchmark:
        ".log/align/filter/{sample}.filter_bam.bm"
    log:
        ".log/align/filter/{sample}.filter_bam.log",
    conda:
        config["conda"]["samtools"]
    # ! 过滤
    # 1.mapq < 20
    # 2.flag 2828:
    #   4(read unmapped)
    #   8(mate unmapped)
    #   256(not primary alignment)
    #   512(read fails platform/vendor quality checks)
    #   2048(supplementary alignment).
    params:
        "-hbS -q 20 -F 2828",
    threads: config["threads"]["medium"]
    shell:
        """
        samtools view -@ {threads} {params} {input} -o {output.bam} 2> {log}
        samtools stat {output.bam} | grep ^SN | cut -f 2- > {output.stat} 2>> {log}
        """


# 探针法需要去重, 引物法不需要去重
rule rmdup_bam:
    input:
        rules.filter_bam.output.bam,
    output:
        sort_name=temp("align/rmdup/{sample}.sorted_by_name.bam"),
        fixmate=temp("align/rmdup/{sample}.fixmate.bam"),
        sort_coodi=temp("align/rmdup/{sample}.sorted_by_coodinate.bam"),
        rmdup="align/rmdup/{sample}.rmdup.bam",
        stats="align/rmdup/{sample}.rmdup_stats.txt",
    benchmark:
        ".log/align/rmdup/{sample}.rmdup_bam.bm"
    log:
        ".log/align/rmdup/{sample}.rmdup_bam.log",
    conda:
        config["conda"]["samtools"]
    params:
        sortn="-n",
        fixmate="-r -c -m",
        markdup="-r -s",
    threads: config["threads"]["low"]
    shell:
        """
        samtools sort -@ 32 {params.sortn} {input} -o {output.sort_name} 2> {log}
        samtools fixmate -@ 32 {params.fixmate} {output.sort_name} {output.fixmate} 2>> {log}
        samtools sort -@ 32 {output.fixmate} -o {output.sort_coodi} 2>> {log}
        samtools markdup -@ 32 {params.markdup} -f {output.stats} {output.sort_coodi} {output.rmdup} 2>> {log}
        """
