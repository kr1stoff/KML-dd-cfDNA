import pandas as pd
import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], 'w')

dfs = []
for bam_stats_file in snakemake.input:
    df = pd.read_csv(bam_stats_file)
    df.insert(0, 'Sample', Path(bam_stats_file).stem.split('.')[0])
    dfs.append(df)
pd.concat(dfs).to_csv(snakemake.output[0], index=False, sep='\t')
