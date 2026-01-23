import pandas as pd

# 读取样本表
samples_df = pd.read_table(
    config["samples_tsv"], header=None, names=["sample", "fq1", "fq2"]
)
samples = samples_df["sample"].tolist()


# 创建软连接
rule create_symlinks:
    input:
        fq1=lambda wildcards: samples_df.loc[
            samples_df["sample"] == wildcards.sample, "fq1"
        ].iloc[0],
        fq2=lambda wildcards: samples_df.loc[
            samples_df["sample"] == wildcards.sample, "fq2"
        ].iloc[0],
    output:
        fq1=".rawdata/{sample}_1.fastq.gz",
        fq2=".rawdata/{sample}_2.fastq.gz",
    log:
        ".log/create_symlinks/{sample}.log",
    benchmark:
        ".log/create_symlinks/{sample}.bm"
    run:
        import os

        # 确保目录存在
        os.makedirs(".rawdata", exist_ok=True)
        # 创建符号链接
        os.symlink(input.fq1, output.fq1)
        os.symlink(input.fq2, output.fq2)
