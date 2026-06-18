# KML-dd-cfDNA 仓库指南

## 项目概述

艾吉泰康 UMI 靶向测序低频变异分析流程。Snakemake 驱动，Rust 工具预处理 UMI。
通过直接对目标位点进行深度 mpileup + 碱基计数计算等位基因频率，最终估算 dd-cfDNA%。
**不使用** fgbio UMI consensus（灵敏度不足），也**不使用** LoFreq/bcftools 进行变异检出。

## 关键命令

```bash
# 运行完整分析流程
snakemake --cores 32 --use-conda \
  --config samples_tsv=$PWD/tests/example.tsv \
  --directory /path/to/results --rerun-incomplete --scheduler greedy

# 编译 Rust UMI 过滤工具（部署前执行）
cargo build --release --manifest-path rust/umi_filter/Cargo.toml
```

## 目录结构

- `Snakefile` — Snakemake 主入口，加载 config.yaml 和 rules/*.smk
- `config.yaml` — 主配置文件（conda 环境映射、线程数、数据库路径）
- `config.schema.yaml` — 配置文件的 JSON Schema 校验
- `rules/` — 7 个规则模块：rawdata, fastq_umi, fastqc, alignment, bam_stats, allele_freq, summary
- `scripts/` — 9 个 Python 脚本（umi_stats/fastp/bam 汇总、碱基解析、AF 计算、dd-cfDNA% 计算、QC 整合）
- `rust/umi_filter/` — Rust UMI 过滤工具，纯 std 实现，无外部依赖
- `assets/` — 资源文件：umi.txt（64 种 UMI 序列）、loci.vcf（500 SNP 位点）、loci.position.txt（samtools 格式）、probeCov.predict.bed（探针覆盖区域）
- `tests/` — 测试数据和脚本
- `desprecated/` — 废弃的旧版流程（fgbio UMI consensus + LoFreq）

## 流程模块

rawdata → fastq_umi → fastqc → alignment → bam_stats → allele_freq → summary

| 模块 | 功能 |
|------|------|
| rawdata | 根据样本 TSV 创建原始 FASTQ 软链接到 `.rawdata/` |
| fastq_umi | pigz 解压 → Rust UMI 过滤（64×64 组合 + Hamming 容错）→ UMI 统计汇总 |
| fastqc | FastQC → MultiQC 聚合 / fastp 质控过滤（UMI 过滤后）→ fastp 统计汇总 |
| alignment | bwa mem（RG 标记）→ sort/index → GATK BQSR（3 个已知位点集）→ sort/index |
| bam_stats | samtools stats/depth（目标区域 + 全基因组）→ 解析为 CSV → BAM 统计汇总 |
| allele_freq | clipOverlap → sort/index → 并行 samtools mpileup（max-depth 500k, min-BQ 30）→ parse_bases 碱基计数 → calc_af 等位基因频率 → calc_ddcfDNA_pct（CareDx 均值法 / GMM） |
| summary | 整合 FASTQ 质控 + BAM 统计 → panel-qc-summary.tsv/.xlsx |

## 架构要点

- **输入**：样本 TSV（无表头，列名：sample, fq1, fq2）
- **参考基因组**：hg19（已在 config.yaml 中配置路径，非 hg38）
- **Rust 工具**接收 6 个参数：`<r1> <r2> <out_r1> <out_r2> <umi_list> <stats>`
  - UMI 为 6bp×2（R1[0:6] + R2[0:6] = 12bp），64 种序列 → 4096 种双端组合
  - 白名单扩增：每个有效 UMI 组合生成所有 Hamming 距离 ≤1 的邻居
  - 通过后截掉序列前 6bp，剩余序列输出
  - 输出统计文件（key: value 格式）：total_reads, kept_reads, total_bases, Q20/Q30 碱基数
- **无需变异检出**：直接对 500 个已知 SNP 位点进行 `samtools mpileup` → `parse_bases.py` 逐位置计数 A/T/C/G → `calc_af.py` 计算 ALT 频率
- **dd-cfDNA% 计算**：`calc_ddcfdfna.py` 支持 CareDx 均值法（默认，乘数 2.11）和 GMM 聚类法
- **`shell.prefix("set +eu;")`**，Snakemake shell 命令前自动添加
- **最终产出**：`qc/multiqc/`（HTML 报告）、`upload/all.dd_cfDNA_pct.xlsx`、`upload/panel-qc-summary.tsv`

## Conda 环境

| 环境名 | 包含工具 |
|--------|---------|
| `basic` | fastqc, multiqc, bwa, samtools, gatk4 |
| `basic2` | fastp, bamutils |
| `python3.12` | python 3.12 + pandas + openpyxl + numpy + scikit-learn + matplotlib + seaborn |

对应 config.yaml 中 `conda:` 下各 key 映射的环境名称，需预先创建。

## 环境依赖

- Snakemake + Conda（各工具在 config.yaml 中按名称映射 conda 环境）
- Rust（仅编译 umi_filter 时需要）
- 参考基因组 hg19 需提前建立 bwa/samtools/picard/GATK 索引
- GATK 已知位点 VCF（1000G、dbSNP、Mills）需提前准备并配置路径
