#!/usr/bin/env python3
from itertools import product
import sys

# 用户配置
# 先解压, python gzip 速度太慢
# pigz -dc SURFS15126-RD-SURF5000-T3113V1-250904-HH05-20251125-ZQJ03_R1.fastq.gz > raw.r1.fq

R1_in = sys.argv[1]
R2_in = sys.argv[2]
R1_out = sys.argv[3]
R2_out = sys.argv[4]
umi_file = sys.argv[5]
umi_stats_outfile = sys.argv[6]

UMI_list = []
with open(umi_file, 'rt') as f:
    for line in f:
        UMI_list.append(line.strip())


MAX_MISMATCH = 1
NUCS = ['A', 'C', 'G', 'T']

# -----------------------------
# 生成 whitelist
# -----------------------------
whitelist = set(u1 + u2 for u1, u2 in product(UMI_list, repeat=2))

# -----------------------------
# 生成所有允许的邻居（Hamming ≤1）
# key = 序列，value = 原始合法 UMI组合
# -----------------------------
allowed = {}
for w in whitelist:
    allowed[w] = w  # exact match
    if MAX_MISMATCH >= 1:
        # Hamming distance = 1
        for i, c in enumerate(w):
            for n in NUCS:
                if n != c:
                    neighbor = w[:i] + n + w[i+1:]
                    allowed[neighbor] = w  # 纠正到原始合法 UMI

print(f"允许序列总数（含 1 mismatch 邻居）: {len(allowed)}")

# -----------------------------
# 处理 FASTQ
# -----------------------------


def count_q(qual_str, qval: int = 20):
    """统计序列中 >= qval 的碱基数, phred33 编码的质量值"""
    return sum((ord(c)-33) >= qval for c in qual_str)


def filter_fastq(R1_in, R2_in, R1_out, R2_out):
    kept_reads, total_reads, total_bases, q20_count, q30_count = 0, 0, 0, 0, 0
    with open(R1_in, 'rt') as f1, open(R2_in, 'rt') as f2, \
            open(R1_out, 'wt') as out1, open(R2_out, 'wt') as out2:

        while True:
            r1_lines = [f1.readline().strip() for _ in range(4)]
            r2_lines = [f2.readline().strip() for _ in range(4)]
            if not r1_lines[0]:
                break
            q20_count += (count_q(r1_lines[3]) + count_q(r2_lines[3]))
            q30_count += (count_q(r1_lines[3], qval=30) +
                          count_q(r2_lines[3], qval=30))
            total_reads += 1
            total_bases += len(r1_lines[1]) + len(r2_lines[1])
            umi1 = r1_lines[1][:6]
            umi2 = r2_lines[1][:6]
            combined = umi1 + umi2

            if combined in allowed:
                out1.write("\n".join(r1_lines) + "\n")
                out2.write("\n".join(r2_lines) + "\n")
                kept_reads += 1

    with open(umi_stats_outfile, 'wt') as f:
        f.write(f"total_reads: {total_reads}\n")
        f.write(f"kept_reads: {kept_reads}\n")
        f.write(f"total_bases: {total_bases}\n")
        f.write(f"q20: {q20_count/total_bases:.4f}\n")
        f.write(f"q30: {q30_count/total_bases:.4f}\n")


# -----------------------------
# 主程序
# -----------------------------
if __name__ == "__main__":
    filter_fastq(R1_in, R2_in, R1_out, R2_out)
