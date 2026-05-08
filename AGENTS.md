# KML-dd-cfDNA 仓库指南

## 项目概述

艾吉泰康 UMI 靶向测序低频变异分析流程。Snakemake 驱动，Rust 工具预处理 UMI。

## 关键命令

```bash
# 运行完整分析流程
snakemake --cores 32 --snakefile wf-igt-500snp/Snakefile --use-conda \
  --config samples_tsv=$PWD/wf-igt-500snp/tests/example.tsv \
  --directory /path/to/results --rerun-incomplete --scheduler greedy

# 编译 Rust UMI 过滤工具（部署前执行）
cd rust/umi_filter && cargo build --release
```

## 目录结构

- `wf-igt-500snp/` — Snakemake 主流程，含 Snakefile、rules/、scripts/、config.yaml
- `rust/umi_filter/` — Rust UMI 过滤工具，独立于 Snakemake 编译
- `desprecated/` — 废弃的旧版流程（UMI consensus + bcftools）

## 架构要点

- 输入：样本 TSV（无表头，列名：sample, fq1, fq2）
- 参考基因组：hg19（已在 config.yaml 中配置路径，非 hg38）
- Rust 工具接收 6 个参数：`<r1> <r2> <out_r1> <out_r2> <umi_list> <stats>`
- 流程模块：rawdata → fastq_umi → fastqc → alignment → bam_stats → variant → allele_freq → summary
- 变异检出使用 **LoFreq**（低频变异专用），而非 bcftools
- 分析策略：用原始数据直接分析，**不使用** fgbio UMI consensus（灵敏度不足）
- UMI 为 6bp，4096 种组合（64×64 双端）
- `shell.prefix("set +eu;")`，Snakemake shell 命令前自动添加

## 环境依赖

- Snakemake + Conda（各工具在 config.yaml 中按名称映射 conda 环境）
- Rust（仅编译 umi_filter 时需要）
- 参考基因组 hg19 需提前建立 bwa/samtools/picard 索引
