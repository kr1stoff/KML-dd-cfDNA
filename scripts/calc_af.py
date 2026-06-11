# calc_af.py
# 功能：根据 loci.anno.vcf 的参考碱基和替换碱基信息，
#       计算 ZQJ03.base_stats 中各 SNP 位点的替换碱基频率（文件1），
#       并利用非 ref/alt 碱基估计背景污染频率（文件2）。
#
# 输入：
#   loci.anno.vcf     — 已知 SNP 位点列表（chr, pos, rsid, ref, alt）
#   ZQJ03.base_stats  — 各染色体位置处 A/T/C/G 碱基计数的统计表
# 输出：
#   ZQJ03.alt_freq.tsv       — 附加了替换碱基频率和背景污染频率列的统计表

import sys

sys.stderr = open(snakemake.log[0], "w")

# 文件路径
vcf_file = snakemake.input.vcf
stats_file = snakemake.input.base_stat
out_file = snakemake.output.alt_freq

# ---------- 第1步：读取 VCF，建立 (染色体, 位置) → (参考碱基, 替换碱基) 的映射 ----------
loci = {}
with open(vcf_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        # VCF 列：chr, pos, rsid, ref, alt （跳过 rsid）
        chrom, pos_str, _, ref, alt = parts[:5]
        loci[(chrom, pos_str)] = (ref.upper(), alt.upper())

# ---------- 第2步：碱基顺序索引 ----------
# base_stats 中 A/T/C/G 四列固定顺序
BASE_ORDER = ['A', 'T', 'C', 'G']
BASE_IDX = {b: i for i, b in enumerate(BASE_ORDER)}

# ---------- 第3步：逐行读取 base_stats，计算并输出 ----------
with open(stats_file) as f, open(out_file, 'w') as o:

    # 读取表头
    header = f.readline().strip()

    # ---- 表头：保留原有列，追加 ALT_BASE, ALT_FREQ, CONTAM_COUNT, CONTAM_FREQ ----
    o.write(header + "\tALT_BASE\tALT_FREQ\tCONTAM_COUNT\tCONTAM_FREQ\n")

    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue

        # 解析 base_stats 的每一行
        chrom, pos, ref = parts[0], parts[1], parts[2].upper()
        counts = [int(parts[3]), int(parts[4]), int(
            parts[5]), int(parts[6])]  # A, T, C, G
        total = int(parts[7])

        key = (chrom, pos)

        if key in loci:
            # ---------- 该位点在 VCF 中已知 ----------
            _, alt = loci[key]
            alt = alt.upper()
            alt_idx = BASE_IDX[alt]
            ref_idx = BASE_IDX[ref]

            # 替换碱基频率
            alt_count = counts[alt_idx]
            alt_freq = alt_count / total if total > 0 else 0

            # 背景污染 = 既不是 ref 也不是 alt 的碱基
            contam_count = sum(counts[i] for i in range(
                4) if i != ref_idx and i != alt_idx)
            contam_freq = contam_count / total if total > 0 else 0

            o.write(
                f"{line.strip()}\t{alt}\t{alt_freq:.6f}\t{contam_count}\t{contam_freq:.6f}\n")

        else:
            # ---------- 该位点不在 VCF 中（理论上不会出现） ----------
            # 全部非 ref 碱基视为背景污染
            ref_idx = BASE_IDX[ref]
            contam_count = total - counts[ref_idx]
            contam_freq = contam_count / total if total > 0 else 0
            o.write(f"{line.strip()}\t.\t.\t{contam_count}\t{contam_freq:.6f}\n")
