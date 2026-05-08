from pathlib import Path
import pandas as pd
import sys

# 将标准错误输出重定向到Snakemake日志文件
sys.stderr = open(snakemake.log[0], "w")

# 初始化字典用于存储所有样本的统计信息
stats_dict = {}

# 遍历所有输入的统计文件
for stats_file in snakemake.input:
    # 从文件名中提取样本名称（去掉.umi.stats后缀）
    sample_name = Path(stats_file).stem.replace('.umi.stats', '')
    stats_dict[sample_name] = {}
    # 读取统计文件内容
    with open(stats_file, 'r') as f:
        for line in f:
            # 按冒号分割键值对
            key, value = line.strip().split(': ')
            stats_dict[sample_name][key] = value

# 将字典转换为DataFrame，以样本名为行索引
df = pd.DataFrame.from_dict(stats_dict, orient='index')
# 设置索引列名
df.index.name = 'Sample'
# 保存为TSV文件
df.to_csv(snakemake.output[0], sep='\t')
