# KML-dd-cfDNA

## 命令行

### 运行

```bash
snakemake --cores 32 --snakefile wf-igt-500snp/Snakefile --use-conda --config samples_tsv=$PWD/tests/example.tsv --directory /data/mengxf/Task/KML260115-dd-cfDNA-iGT/results/260130 --rerun-incomplete --scheduler greedy
```

## 开发说明

20260228

- 不再使用 fgbio umi consensus 方式, 使用原始数据直接分析
  fgbio umi consensus 分析流程灵敏度无法达到 < 1%, 1% VAF 以下非常不准确
- 不再使用 `bcftools` 算 VCF, 改用低频变异专用的 `Mutect2`
  IGV 核对原始比对 BAM 时发现, 很多 10% VAF 的变异都被

20260127

- 艾吉泰康UMI
  IGT_UMI_Adapter_and_UDI_Primer 是 64*64=4096 种 UMI 组合, UMI 长度为 6bp
