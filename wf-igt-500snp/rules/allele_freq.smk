rule mpileup_bcftools:
    input:
        bam=rules.sort_umi_bam.output.bam,
        ref=config["database"]["hg19"],
        region=config["database"]["region_slop500bp"],
    output:
        vcf="snp/umi/{sample}.region.call.vcf.gz",
        csi="snp/umi/{sample}.region.call.vcf.gz.csi",
    benchmark:
        ".log/snp/umi/{sample}.mpileup_bcftools.bm"
    log:
        ".log/snp/umi/{sample}.mpileup_bcftools.log",
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


rule extract_500snps:
    input:
        vcf=rules.mpileup_bcftools.output.vcf,
        anno=config["database"]["anno_500snp"],
    output:
        txt="snp/umi/{sample}.500snps.txt",
    benchmark:
        ".log/snp/umi/{sample}.extract_500snps.bm"
    log:
        ".log/snp/umi/{sample}.extract_500snps.log",
    conda:
        config["conda"]["bcftools"]
    threads: config["threads"]["medium"]
    shell:
        """
        bcftools view -R {input.anno} {input.vcf} 2> {log} | \
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%AD]\n' > {output.txt} 2>> {log}
        """
