#!/usr/bin/env python3
from itertools import product
import sys

# 用户配置
# 先解压, python gzip 速度太慢
# pigz -dc SURFS15126-RD-SURF5000-T3113V1-250904-HH05-20251125-ZQJ03_R1.fastq.gz > raw.r1.fq

# R1_in = "raw.r1.fq"
# R2_in = "raw.r2.fq"
R1_in = sys.argv[1]
R2_in = sys.argv[2]
R1_out = "filtered.r1.fq"
R2_out = "filtered.r2.fq"

UMI_list = [
    'TCGGTA', 'GACGTT', 'AAGTTG', 'AAGGCA', 'AAGCGT', 'AACTGA', 'AACCAG', 'ATAAGC', 'ATACCT', 'ATTGAG',
    'ATGCAA', 'ATCATG', 'ATCTAC', 'AGATGT', 'AGACTC', 'AGTACA', 'AGTGGC', 'ACAACG', 'ACAGAT', 'TAATCG',
    'TATGAC', 'TATCCA', 'TACTTC', 'TTACGA', 'TTGTCA', 'TTCACT', 'TTCTGG', 'TGAGCA', 'TGTATC', 'TGTTAG',
    'TGTCGT', 'TGCCAC', 'TCAGGC', 'TCTTGA', 'TCTGCT', 'TCGACC', 'TCGCAG', 'TCCTAT', 'GAATAC', 'GATGGA',
    'GACAGC', 'GTAGAA', 'GTTAAC', 'GTTCTA', 'GTGAGT', 'GGAACC', 'GGATTA', 'GGCAAG', 'GCAATT', 'GCACCA',
    'GCTTCG', 'GCTCAT', 'GCGTAA', 'CAAGTA', 'CAACAT', 'CATTGT', 'CATGCG', 'CTAGGT', 'CTGATA', 'CTGGAC',
    'CTCCTC', 'CGTCTG', 'CCATAG', 'CCGCTT'
]

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


def filter_fastq(R1_in, R2_in, R1_out, R2_out):
    kept, total = 0, 0
    with open(R1_in, 'rt') as f1, open(R2_in, 'rt') as f2, \
            open(R1_out, 'wt') as out1, open(R2_out, 'wt') as out2:

        while True:
            r1_lines = [f1.readline().strip() for _ in range(4)]
            r2_lines = [f2.readline().strip() for _ in range(4)]
            if not r1_lines[0]:
                break

            total += 1
            umi1 = r1_lines[1][:6]
            umi2 = r2_lines[1][:6]
            combined = umi1 + umi2

            if combined in allowed:
                out1.write("\n".join(r1_lines) + "\n")
                out2.write("\n".join(r2_lines) + "\n")
                kept += 1

    print(f"总 reads: {total}, 保留 reads: {kept}, 过滤掉: {total-kept}")


# -----------------------------
# 主程序
# -----------------------------
if __name__ == "__main__":
    filter_fastq(R1_in, R2_in, R1_out, R2_out)
