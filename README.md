# KML-dd-cfDNA

## 命令行

### 运行

- 主流程

  ```bash
  snakemake --cores 32 --use-conda --rerun-incomplete --scheduler greedy \
    --config samples_tsv=$PWD/wf-igt-500snp/tests/example.tsv \
    --directory /data/mengxf/Task/KML260115-dd-cfDNA-iGT/results/260130
  ```

s

- 编译 Rust 程序 (分析流程部署前)
编译后可执行文件软连接到 `wf-igt-500snp\tools`

  ```bash
  cd rust/umi_filter
  cargo build --release
  ```

## 开发说明
- 20260509
  - 不再使用 GATK MarkDuplicates 步骤, 重复序列被标记为 1024(0x400) 下后 LoFreq 和 IGV 都会默认忽略掉重复序列.
  - 比对后使用 BQSR 校正

- 20260317
  - 使用 Rust 处理原始数据, 过滤不含 UMI 的读取, 并截掉 UMI 部分

- 20260228
  - 不再使用 fgbio umi consensus 方式, 使用原始数据直接分析, fgbio umi consensus 分析流程灵敏度无法达到 < 1%, 1% VAF 以下非常不准确
  - 不再使用 `bcftools` 算 VCF, 改用 `LoFreq`
  - IGV 核对原始比对 BAM 时发现, 很多 10% VAF 的变异都被过滤掉了

- 20260127
  - 艾吉泰康UMI, IGT_UMI_Adapter_and_UDI_Pimer 是 64*64=4096 种 UMI 组合, UMI 长度为 6bp
