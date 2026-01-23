rule bwa_mem_hg19:
    input:
        rules.fastp.output.trimmed,
    output:
        temp("align/bwa/{sample}.bam"),
    benchmark:
        ".log/align/bwa/{sample}.hg19.bm"
    log:
        ".log/align/bwa/{sample}.hg19.log",
    conda:
        config["conda"]["bwa"]
    params:
        bwa="-M -Y -R '@RG\\tID:{sample}\\tSM:{sample}'",
        view="-hbS",
    threads: config["threads"]["medium"]
    shell:
        """
        bwa mem -t {threads} {params.bwa} {config[database][hg19]} {input} 2> {log} | \
            samtools view -@ {threads} {params.view} - 2>> {log} | \
            samtools sort -@ {threads} -o {output} - 2>> {log}
        """
