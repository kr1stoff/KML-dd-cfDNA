import pandas as pd
import sys

# 将标准错误输出重定向到Snakemake日志文件
sys.stderr = open(snakemake.log[0], "w")

umi_stats_file = snakemake.input.umi[0]
fastp_stats_file = snakemake.input.fastp[0]

df_umi = pd.read_csv(umi_stats_file, sep='\t', index_col='Sample')
df_fastp = pd.read_csv(fastp_stats_file, sep='\t', index_col='Sample')

df_merged = pd.merge(df_umi, df_fastp, left_index=True, right_index=True)
df_merged.to_csv(snakemake.output[0], sep='\t')
