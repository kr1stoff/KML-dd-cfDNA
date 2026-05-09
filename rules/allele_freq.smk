rule bcftools_view_500snps:
    input:
        vcf=rules.variant_bgzip_and_index.output.vcf,
        anno=f"{workflow.basedir}/assets/loci.anno.vcf",
    output:
        txt="variant/snp/{sample}.500snps.raw.txt",
    benchmark:
        ".log/variant/snp/{sample}.bcftools_view_500snps.bm"
    log:
        ".log/variant/snp/{sample}.bcftools_view_500snps.log",
    conda:
        config["conda"]["bcftools"]
    threads: config["threads"]["medium"]
    shell:
        """
        bcftools view -R {input.anno} {input.vcf} 2> {log} | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AF\n' > {output.txt} 2>> {log}
        """


rule extract_500snps:
    input:
        txt=rules.bcftools_view_500snps.output.txt,
        anno=f"{workflow.basedir}/assets/loci.anno.vcf",
    output:
        txt="variant/snp/{sample}.500snps.txt",
    benchmark:
        ".log/variant/snp/{sample}.extract_500snps.bm"
    log:
        ".log/variant/snp/{sample}.extract_500snps.log",
    conda:
        config["conda"]["python"]
    run:
        targets = []
        with open(input.anno, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                del fields[2]  # 删除第3个元素
                targets.append(fields)
        with open(input.txt, "r") as f, open(output.txt, "w") as o:
            for line in f:
                if line.split("\t")[:4] in targets:
                    o.write(line)
