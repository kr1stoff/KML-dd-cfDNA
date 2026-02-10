# KML-dd-cfDNA

## 命令行

### Snakemake 运行

```bash
snakemake --cores 32 --snakefile wf-igt-500snp/Snakefile --use-conda --config samples_tsv=$PWD/tests/example.tsv --directory /data/mengxf/Task/KML260115-dd-cfDNA-iGT/results/260130 --rerun-incomplete
```

## 说明

1. 艾吉泰康UMI
`IGT_UMI_Adapter_and_UDI_Primer` 是 64*64=4096 种 UMI 组合, UMI 长度为 6bp
