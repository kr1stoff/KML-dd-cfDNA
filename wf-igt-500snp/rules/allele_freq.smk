rule extract_500snps_raw:
    input:
        vcf=rules.variant_bgzip_and_index.output.vcf,
        anno=f"{workflow.basedir}/assets/loci.anno.vcf",
    output:
        txt="variant/snp/{sample}.500snps.txt",
    benchmark:
        ".log/variant/snp/{sample}.extract_500snps_raw.bm"
    log:
        ".log/variant/snp/{sample}.extract_500snps_raw.log",
    conda:
        config["conda"]["bcftools"]
    threads: config["threads"]["medium"]
    shell:
        # 总深度不要用 DP, DP 是 raw depth, 不是高质量 depth. 总深度使用 %AD 加和
        """
        bcftools view -R {input.anno} {input.vcf} 2> {log} | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP][\t%AF]\n' > {output.txt} 2>> {log}
        """
