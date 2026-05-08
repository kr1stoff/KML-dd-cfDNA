import pandas as pd
from pathlib import Path
import sys

sys.stderr = open(snakemake.log[0], 'w')


def get_bam_stat_dict(bam_stat_file) -> dict:
    bam_stat_dict = {}
    with open(bam_stat_file, 'r') as f:
        for line in f:
            if line.startswith('SN'):
                lns = line.split('#')[0].strip().split('\t')
                bam_stat_dict[lns[1].replace(':', '')] = lns[2]
    return bam_stat_dict


# 靶区域 bam stat 字典
target_stats = get_bam_stat_dict(snakemake.input[1])
# 全部比对结果 bam stat 字典
all_stats = get_bam_stat_dict(snakemake.input[0])
# 比对率 reads mapped and paired / raw total sequences
mapped_reads = int(all_stats['reads mapped and paired'])
raw_total_sequences = int(all_stats['raw total sequences'])
mapped_rate = mapped_reads / raw_total_sequences if raw_total_sequences > 0 else 0
# 在靶率 (target) reads mapped and paired / (all) reads mapped and paired
ontarget_reads = int(target_stats['reads mapped and paired'])
ontarget_rate = ontarget_reads / mapped_reads if mapped_reads > 0 else 0

# 覆盖率
# ! 阴控或其他没有比对到参考的样本会报错
target_depth_file = Path(snakemake.input[2])
if target_depth_file.stat().st_size != 0:
    target_depth = pd.read_table(target_depth_file, header=None)
    target_size = target_depth.shape[0]
    target_covered_size = target_depth[target_depth[2] > 0].shape[0]
    covarage_rate = target_covered_size / target_size
    # 4X/10X/30X/50X/100X 覆盖深度
    cover4x_rate = target_depth[target_depth[2] > 4].shape[0] / target_size
    cover10x_rate = target_depth[target_depth[2] > 10].shape[0] / target_size
    cover30x_rate = target_depth[target_depth[2] > 30].shape[0] / target_size
    cover100x_rate = target_depth[target_depth[2] > 100].shape[0] / target_size
    cover200x_rate = target_depth[target_depth[2] > 200].shape[0] / target_size
    # 20%/50% 均一性
    mean_depth = int(target_depth[2].mean())
    depth_20_rate = target_depth[target_depth[2] > mean_depth * 0.2].shape[0] / target_size
    depth_50_rate = target_depth[target_depth[2] > mean_depth * 0.5].shape[0] / target_size
    # [20251013] 新增均一性 P90/P10 数值
    depth_p90 = int(target_depth[2].quantile(0.9))
    depth_p10 = max(int(target_depth[2].quantile(0.1)), 1)  # 避免除0错误
    p90_div_p10 = depth_p90 / depth_p10
else:
    target_size, target_covered_size, covarage_rate = 0, 0, 0
    cover4x_rate, cover10x_rate, cover30x_rate, cover100x_rate, cover200x_rate = 0, 0, 0, 0, 0
    mean_depth, depth_20_rate, depth_50_rate, p90_div_p10 = 0, 0, 0, 0

# 输出
outputs = [mapped_reads, mapped_rate, ontarget_reads, ontarget_rate, target_size, target_covered_size, covarage_rate,
           cover4x_rate, cover10x_rate, cover30x_rate, cover100x_rate, cover200x_rate,
           mean_depth, depth_20_rate, depth_50_rate, p90_div_p10]
outputs = [o if type(o) == int else round(o, 6) for o in outputs]
df = pd.DataFrame([outputs], columns=[
    'MappedReads', 'MappingRate', 'OnTargetReads', 'OnTargetRate', 'TargetSize', 'TargetCoveredSize', 'CoverageRate',
    '4xCoverageRate', '10xCoverageRate', '30xCoverageRate', '100xCoverageRate', '200xCoverageRate',
    'MeanDepth', '20%MeanDepthCoverageRate', '50%MeanDepthCoverageRate', 'P90/P10'])
df.to_csv(snakemake.output[0], index=False)
