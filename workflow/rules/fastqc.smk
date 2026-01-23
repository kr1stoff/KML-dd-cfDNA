# QC
rule fastqc:
    input:
        rules.create_symlinks.output.fq1,
        rules.create_symlinks.output.fq2,
    output:
        directory("qc/fastqc/{sample}"),
    benchmark:
        ".log/qc/fastqc/{sample}.bm"
    log:
        ".log/qc/fastqc/{sample}.log",
    conda:
        config["conda"]["fastqc"]
    threads: config["threads"]["low"]
    shell:
        "mkdir {output} && fastqc {input} -o {output} -t {threads} --extract &> {log}"


rule multiqc:
    input:
        expand("qc/fastqc/{sample}", sample=samples),
    output:
        directory("multiqc"),
    benchmark:
        ".log/qc/multiqc/multiqc.bm"
    log:
        ".log/qc/multiqc/multiqc.log",
    conda:
        config["conda"]["multiqc"]
    shell:
        "multiqc {input} --outdir {output} 2> {log}"
