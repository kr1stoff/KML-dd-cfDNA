rule mpileup_bcftools_raw:
    input:
        bam=rules.bwa_mem_raw.output.bam,
        ref=config["database"]["hg19"],
        region=config["database"]["region_slop500bp"],
    output:
        vcf="snp/raw/{sample}.region.call.vcf.gz",
        csi="snp/raw/{sample}.region.call.vcf.gz.csi",
    benchmark:
        ".log/snp/raw/{sample}.mpileup_bcftools_raw.bm"
    log:
        ".log/snp/raw/{sample}.mpileup_bcftools_raw.log",
    conda:
        config["conda"]["bcftools"]
    threads: config["threads"]["medium"]
    params:
        mpileup="--max-depth 25000 --min-MQ 20 --min-BQ 30 --no-BAQ -Ou",
        call="--multiallelic-caller -Ov",
    shell:
        """
        bcftools mpileup --threads {threads} {params.mpileup} --fasta-ref {input.ref} --regions-file {input.region} {input.bam} 2> {log} | \
            bcftools call --threads {threads} {params.call} 2>> {log} | \
            bgzip -c > {output.vcf} 2>> {log}
        bcftools index {output.vcf} 2>> {log}
        """


rule extract_500snps_raw:
    input:
        vcf=rules.mpileup_bcftools_raw.output.vcf,
        anno=config["database"]["anno_500snp"],
    output:
        txt="snp/raw/{sample}.500snps.txt",
    benchmark:
        ".log/snp/raw/{sample}.extract_500snps_raw.bm"
    log:
        ".log/snp/raw/{sample}.extract_500snps_raw.log",
    conda:
        config["conda"]["bcftools"]
    threads: config["threads"]["medium"]
    shell:
        """
        bcftools view -R {input.anno} {input.vcf} 2> {log} | \
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%AD]\n' > {output.txt} 2>> {log}
        """
