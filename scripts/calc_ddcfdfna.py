#!/usr/bin/env python3
"""
dd-cfDNA 占比计算脚本

【算法原理说明】
本脚本实现了两种 dd-cfDNA(供体来源游离 DNA)占比估算方法：

1. CareDx AlloSure 均值法(默认方法):
   - 核心思想：无需供/受者先验基因分型，通过统计推断估算供体 DNA 比例
   - 算法流程:
     1. 过滤低深度位点(TOTAL < 1000)和高污染位点(CONTAM_FREQ > 0.01)
     2. 对每个 SNP，计算 minor_AF = min(ALT_FREQ, 1 - ALT_FREQ)，即次要等位基因频率
     3. 推断受者基因型：minor_AF <= het_threshold → 受者纯合(信息位点)
        - 原理：若受者为 AA 纯合，检测到的 ALT 等位基因只能来自供体
     4. 收集所有信息位点的 donor_AF = minor_AF
     5. 剔除 top 5% 离群值(消除 mapping 伪影、PCR 错误等极端值)
     6. 取均值得到 mean_donor_AF
     7. 乘以亲缘关系校正乘数 M：无关供受者 M=2.11，亲子 M≈1.7，同胞 M≈1.4

2. GMM(高斯混合模型)法(备选方法):
   - 核心思想：通过聚类分析识别等位基因频率分布中的特征峰
   - 适用于较高 dd-cfDNA 比例场景(如 >10%)
   - 模型假设：纯合位点的 minor_AF 分布包含多个高斯成分
     - background: 测序噪声背景(~0)
     - AB: 供体杂合位点贡献(频率≈f/2)
     - BB: 供体纯合位点贡献(频率≈f)
   - 从 BB 聚类中心直接估算 f

【计算方法选择说明】
--method 参数支持三种选择:
  * caredx: 使用 CareDx 均值法(推荐，适用于低比例场景，如 0.1%-10%)
  * gmm: 使用 GMM 聚类法(适用于高比例场景，如 >10%)
  * both: 同时运行两种方法，便于结果对比

【输入输出】
输入: data/ 目录下的 *.alt_freq.tsv 文件(单样本即可计算)
输出: dd-cfDNA% 占比数值及详细质控统计

【使用示例】
    python calc_ddcfdfna.py data/ZQJ09.alt_freq.tsv
    python calc_ddcfdfna.py data/*.alt_freq.tsv --output summary.tsv
    python calc_ddcfdfna.py data/*.alt_freq.tsv --method both

【实验背景】
    500 SNP 探针捕获 + 固定 64x64 UMI(共 4096 种组合)
    测序流程: fastp → bwa mem (hg19) → BQSR → clipOverlap → mpileup → 自编等位基因计数
"""

# import argparse
# import glob
# import warnings
import os
import numpy as np
import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


def compute_ddcfdfna_caredx(tsv_path, multiplier=2.11, min_depth=1000,
                            het_threshold=0.35, outlier_fraction=0.05,
                            max_contam=0.01):
    """CareDx AlloSure 风格的均值法 dd-cfDNA 估算

    【算法原理】
    假设受者为 AA 纯合时，检测到的 ALT 等位基因(如 B)只能来自供体。
    对于每个 SNP，minor_AF <= het_threshold 表明受者很可能是纯合的
    (因为杂合子的两个等位基因频率应接近 0.5)。
    收集所有信息位点的 minor_AF，取均值并乘以亲缘关系校正系数得到最终估算值。

    参数
    ----------
    tsv_path : str
        ALT 频率 TSV 文件路径
    multiplier : float
        亲缘关系校正乘数(无关供受者=2.11，亲子≈1.7，同胞≈1.4)
        - 理论依据：HLA 不匹配时，供体特异性等位基因在受者样本中的最大可能频率
    min_depth : int
        保留位点的最小总测序深度(默认 1000x，确保统计可靠性)
    het_threshold : float
        杂合判断阈值，minor_AF <= 此值认为受者纯合(信息位点)
    outlier_fraction : float
        剔除的顶部离群值比例(默认 5%，消除 mapping 伪影等极端值)
    max_contam : float
        最大可接受的污染频率，超过此值的位点被排除

    返回
    -------
    dict 包含以下键:
        sample: 样本名称
        n_total_raw: 原始位点总数
        n_total: 深度过滤后位点数量
        n_after_contam: 污染过滤后位点数量
        n_homozygous: 推断为受者纯合的信息位点数量
        n_used: 去除离群值后最终使用的位点数量
        mean_donor_AF: 供体等位基因频率均值
        dd_cfDNA_pct: dd-cfDNA 占比(%)
        q5/25/median/75/95_donor_AF: donor_AF 分布的分位数统计
        warning: 质控警告信息
    """
    sample_name = os.path.splitext(os.path.basename(tsv_path))[0]

    # ========== 1. 读取输入文件 ==========
    # 输入文件为 TSV 格式，包含列: CHROM, POS, TOTAL, ALT_FREQ, CONTAM_FREQ
    df = pd.read_csv(tsv_path, sep='\t', comment='#')
    required = ['CHROM', 'POS', 'TOTAL', 'ALT_FREQ', 'CONTAM_FREQ']
    for col in required:
        if col not in df.columns:
            return {'sample': sample_name, 'error': f'Missing column: {col}'}

    n_total_raw = len(df)  # 原始位点总数

    # ========== 2. 质量过滤 ==========
    # 过滤低深度位点(测序深度不足会导致频率估算不可靠)
    df = df[df['TOTAL'] >= min_depth].copy()
    n_total = len(df)

    # 过滤高污染位点(CONTAM_FREQ 指示样本污染程度)
    df = df[df['CONTAM_FREQ'] < max_contam].copy()
    n_after_contam = len(df)

    # 检查过滤后位点数量是否足够
    if n_after_contam < 20:
        return {'sample': sample_name, 'n_total_raw': n_total_raw,
                'n_total': n_total, 'n_after_contam': n_after_contam,
                'error': f'Too few sites after filtering ({n_after_contam})'}

    # ========== 3. 计算次要等位基因频率 ==========
    # minor_AF = min(ALT_FREQ, 1-ALT_FREQ)
    # 将等位基因频率统一转换到 [0, 0.5] 区间
    minor_AF = np.where(df['ALT_FREQ'].values <= 0.5,
                        df['ALT_FREQ'].values,
                        1.0 - df['ALT_FREQ'].values)

    # ========== 4. 推断受者纯合位点 ==========
    # 假设：minor_AF <= het_threshold → 受者纯合(信息位点)
    # 原理：杂合子的两个等位基因频率应接近 0.5，纯合子的次要等位基因频率应很低
    hom_mask = minor_AF <= het_threshold
    n_homozygous = int(hom_mask.sum())

    # 检查纯合位点数是否足够
    if n_homozygous < 10:
        return {'sample': sample_name, 'n_total_raw': n_total_raw,
                'n_total': n_total, 'n_after_contam': n_after_contam,
                'n_homozygous': n_homozygous,
                'error': f'Too few homozygous sites ({n_homozygous})'}

    # 收集信息位点的供体等位基因频率
    donor_AF = minor_AF[hom_mask]

    # ========== 5. 离群值去除 ==========
    # 剔除顶部 5% 的极端值，消除 mapping 伪影、PCR 错误等影响
    if outlier_fraction > 0 and len(donor_AF) > 0:
        cutoff = np.quantile(donor_AF, 1.0 - outlier_fraction)
        donor_AF = donor_AF[donor_AF <= cutoff]
    n_used = len(donor_AF)

    if n_used < 5:
        return {'sample': sample_name, 'n_total_raw': n_total_raw,
                'n_total': n_total, 'n_after_contam': n_after_contam,
                'n_homozygous': n_homozygous, 'n_used': n_used,
                'error': 'Too few sites after outlier removal'}

    # ========== 6. 计算最终结果 ==========
    # 取均值作为供体等位基因频率估算
    mean_donor_AF = float(np.mean(donor_AF))
    # 乘以亲缘关系校正乘数得到 dd-cfDNA 占比(转换为百分比)
    dd_cfDNA_pct = mean_donor_AF * multiplier * 100.0

    # ========== 7. 质控警告 ==========
    warnings = []
    if n_homozygous < 60:
        warnings.append(f'low_homozygous_count({n_homozygous})')  # 信息位点过少
    if dd_cfDNA_pct < 0.05:
        warnings.append('below_typical_LOD(0.1%)')  # 低于典型检测下限
    if dd_cfDNA_pct > 20:
        warnings.append('above_model_limit')  # 超出模型适用范围

    return {
        'sample': sample_name,
        'method': 'caredx_mean',
        'n_total_raw': n_total_raw,
        'n_total': n_total,
        'n_after_contam': n_after_contam,
        'n_homozygous': n_homozygous,
        'n_used': n_used,
        'mean_donor_AF': mean_donor_AF,
        'dd_cfDNA_pct': dd_cfDNA_pct,
        'q5_donor_AF': float(np.quantile(donor_AF, 0.05)),
        'q25_donor_AF': float(np.quantile(donor_AF, 0.25)),
        'median_donor_AF': float(np.median(donor_AF)),
        'q75_donor_AF': float(np.quantile(donor_AF, 0.75)),
        'q95_donor_AF': float(np.quantile(donor_AF, 0.95)),
        'warning': '; '.join(warnings),
    }


def compute_ddcfdfna_gmm(tsv_path, het_threshold=0.35, max_contam=0.01,
                         min_depth=1000):
    """GMM(高斯混合模型)法 dd-cfDNA 估算(适用于较高比例场景)

    【算法原理】
    通过高斯混合模型对纯合位点的 minor_AF 分布进行聚类分析，识别特征峰：
    - background: 测序噪声背景，频率接近 0
    - AB: 供体杂合位点贡献，频率≈f/2(供体提供一个等位基因)
    - BB: 供体纯合位点贡献，频率≈f(供体提供两个等位基因)

    从 BB 聚类中心直接估算 dd-cfDNA 比例 f，无需亲缘关系校正乘数。

    【适用场景】
    当 dd-cfDNA 比例较高(如 >10%)时，均值法可能受离群值影响较大，
    GMM 法通过聚类可以更准确地识别真实的供体信号。

    参数
    ----------
    tsv_path : str
        ALT 频率 TSV 文件路径
    het_threshold : float
        杂合判断阈值
    max_contam : float
        最大污染频率过滤
    min_depth : int
        最小测序深度

    返回
    -------
    dict 包含聚类结果和 dd-cfDNA 估算值
    """
    try:
        from sklearn.mixture import GaussianMixture
    except ImportError:
        return {'error': 'sklearn not available; install: pip install scikit-learn'}

    sample_name = os.path.splitext(os.path.basename(tsv_path))[0]

    df = pd.read_csv(tsv_path, sep='\t', comment='#')
    df = df[(df['TOTAL'] >= min_depth) & (
        df['CONTAM_FREQ'] < max_contam)].copy()

    if len(df) < 20:
        return {'sample': sample_name, 'error': 'Too few sites after filtering'}

    minor_AF = np.where(df['ALT_FREQ'].values <= 0.5,
                        df['ALT_FREQ'].values,
                        1.0 - df['ALT_FREQ'].values)

    hom_mask = minor_AF <= het_threshold
    donor_AF = minor_AF[hom_mask].reshape(-1, 1)

    if len(donor_AF) < 20:
        return {'sample': sample_name, 'n_homozygous': int(hom_mask.sum()),
                'error': 'Too few homozygous sites for GMM'}

    # ========== GMM 聚类分析 ==========
    # 尝试 2-4 个聚类成分，通过 BIC(贝叶斯信息准则)选择最佳模型
    # BIC 越小表示模型拟合优度越高(同时惩罚过多参数)
    best_bic = np.inf
    best_gmm = None
    best_k = 2
    for k in range(2, 5):
        gmm = GaussianMixture(n_components=k, random_state=0, n_init=5)
        gmm.fit(donor_AF)
        bic = gmm.bic(donor_AF)
        if bic < best_bic:
            best_bic = bic
            best_gmm = gmm
            best_k = k

    # 获取排序后的聚类中心均值
    means = np.sort(best_gmm.means_.flatten())

    # ========== 从聚类结果估算 dd-cfDNA 比例 ==========
    # 模型假设：聚类中心对应 [background, AB, BB] 三个成分
    # background: 测序噪声(频率<0.001)
    # AB: 供体杂合位点贡献(频率≈f/2)
    # BB: 供体纯合位点贡献(频率≈f)

    # 过滤掉背景噪声成分
    non_bg = means[means > 0.001]

    if len(non_bg) >= 2:
        # 有 AB 和 BB 两个信号峰
        f_est = non_bg[-1]      # BB 峰 → f
        f_half_est = non_bg[-2]  # AB 峰 → f/2
    elif len(non_bg) == 1:
        # 只有一个信号峰(可能是 BB 或 AB)
        f_est = non_bg[0]
        f_half_est = None
    else:
        # 没有有效信号
        f_est = 0.0
        f_half_est = None

    return {
        'sample': sample_name,
        'method': 'gmm',
        'n_homozygous': int(hom_mask.sum()),
        'n_clusters': best_k,
        'cluster_means': [round(float(m), 6) for m in means],
        'f_estimate_gmm': float(f_est),
        'f_half_estimate': float(f_half_est) if f_half_est else None,
        'dd_cfDNA_pct_gmm': float(f_est * 100.0),
        'dd_cfDNA_pct_from_AB': float(f_half_est * 2 * 100.0) if f_half_est else None,
    }


############################################ Main ############################################
# IO
input_files = snakemake.input
output_tsv = snakemake.output.tsv
output_xlsx = snakemake.output.xlsx

# Parameters
method = snakemake.params.method
multiplier = snakemake.params.multiplier
min_depth = snakemake.params.min_depth
het_threshold = snakemake.params.het_threshold
outlier_fraction = snakemake.params.outlier_fraction
max_contam = snakemake.params.max_contam

# Run calculation
results = []

for fpath in input_files:
    if method in ('caredx', 'both'):
        r = compute_ddcfdfna_caredx(
            fpath,
            multiplier=multiplier,
            min_depth=min_depth,
            het_threshold=het_threshold,
            outlier_fraction=outlier_fraction,
            max_contam=max_contam,
        )
        results.append(r)
        if 'error' in r:
            print(
                f'{os.path.basename(fpath):20s}  ERROR: {r["error"]}', file=sys.stderr)
        else:
            print(f'{r["sample"]:20s}  dd-cfDNA = {r["dd_cfDNA_pct"]:8.4f}%  '
                  f'(hom={r["n_homozygous"]:>3d}, used={r["n_used"]:>3d}, '
                  f'mean_AF={r["mean_donor_AF"]:.6f})'
                  + (f'  [{r["warning"]}]' if r['warning'] else ''), file=sys.stderr)

    if method in ('gmm', 'both'):
        r2 = compute_ddcfdfna_gmm(
            fpath,
            het_threshold=het_threshold,
            max_contam=max_contam,
            min_depth=min_depth,
        )
        if 'error' not in r2:
            results.append(r2)
            print(f'{r2["sample"]:20s}  [GMM] dd-cfDNA = {r2["dd_cfDNA_pct_gmm"]:8.4f}%  '
                  f'(clusters={r2["cluster_means"]})', file=sys.stderr)
        else:
            print(
                f'{os.path.basename(fpath):20s}  [GMM] ERROR: {r2["error"]}', file=sys.stderr)

if results:
    df_out = pd.DataFrame(results).sort_values(by='sample')
    df_out.to_csv(output_tsv, sep='\t', index=False)
    df_out.to_excel(output_xlsx, index=False)
    print(f'\nSummary saved to: {output_tsv}', file=sys.stderr)
else:
    # 没有有效结果，输出空文件
    print(
        f'[WARNING] No valid results to save to {output_tsv} and {output_xlsx}', file=sys.stderr)
    pd.DataFrame().to_csv(output_tsv, sep='\t', index=False)
    pd.DataFrame().to_excel(output_xlsx, index=False)

"""
def main():
    parser = argparse.ArgumentParser(
        description='dd-cfDNA% 计算 — 从单个时间点 ALT 频率 TSV 估算供体来源 DNA 占比')
    parser.add_argument('input', nargs='+',
                        help='TSV file(s) or glob pattern')
    parser.add_argument('--method', choices=['caredx', 'gmm', 'both'],
                        default='caredx',
                        help='计算方法选择: caredx: CareDx 均值法(推荐，适用于低比例 0.1%%-10%%); gmm: GMM 聚类法(适用于高比例 >10%%); both: 同时运行两种方法进行对比')
    parser.add_argument('--multiplier', type=float, default=2.11,
                        help='亲缘关系校正乘数 M (default: 2.11, 无关供受者)')
    parser.add_argument('--min-depth', type=int, default=1000,
                        help='最小总深度 (default: 1000)')
    parser.add_argument('--het-threshold', type=float, default=0.35,
                        help='杂合判断阈值 (default: 0.35)')
    parser.add_argument('--outlier-fraction', type=float, default=0.05,
                        help='剔除 top 离群值比例 (default: 0.05)')
    parser.add_argument('--max-contam', type=float, default=0.01,
                        help='最大 CONTAM_FREQ 过滤 (default: 0.01)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='输出汇总 TSV 文件路径')

    args = parser.parse_args()

    results = []

    for pattern in args.input:
        matched = sorted(glob.glob(pattern))
        if not matched:
            print(f'[WARNING] No files match: {pattern}', file=sys.stderr)
            continue
        for fpath in matched:
            if not os.path.isfile(fpath):
                continue

            if args.method in ('caredx', 'both'):
                r = compute_ddcfdfna_caredx(
                    fpath,
                    multiplier=args.multiplier,
                    min_depth=args.min_depth,
                    het_threshold=args.het_threshold,
                    outlier_fraction=args.outlier_fraction,
                    max_contam=args.max_contam,
                )
                results.append(r)
                if 'error' in r:
                    print(f'{os.path.basename(fpath):20s}  ERROR: {r["error"]}')
                else:
                    print(f'{r["sample"]:20s}  dd-cfDNA = {r["dd_cfDNA_pct"]:8.4f}%  '
                          f'(hom={r["n_homozygous"]:>3d}, used={r["n_used"]:>3d}, '
                          f'mean_AF={r["mean_donor_AF"]:.6f})'
                          + (f'  [{r["warning"]}]' if r['warning'] else ''))

            if args.method in ('gmm', 'both'):
                r2 = compute_ddcfdfna_gmm(
                    fpath,
                    het_threshold=args.het_threshold,
                    max_contam=args.max_contam,
                    min_depth=args.min_depth,
                )
                if 'error' not in r2:
                    results.append(r2)
                    print(f'{r2["sample"]:20s}  [GMM] dd-cfDNA = {r2["dd_cfDNA_pct_gmm"]:8.4f}%  '
                          f'(clusters={r2["cluster_means"]})')
                else:
                    print(f'{os.path.basename(fpath):20s}  [GMM] ERROR: {r2["error"]}')

    if args.output and results:
        df_out = pd.DataFrame(results)
        df_out.to_csv(args.output, sep='\t', index=False)
        print(f'\nSummary saved to: {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
"""
