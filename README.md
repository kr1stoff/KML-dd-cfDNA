# KML-dd-cfDNA

## 命令行

### Snakemake 运行

```bash
snakemake --snakefile wf-igt-500snp/Snakefile --cores 32 --use-conda --config samples_tsv=/data/mengxf/Task/KML260115-dd-cfDNA-iGT/work/260123-input/input.tsv --directory /data/mengxf/Task/KML260115-dd-cfDNA-iGT/results/250126 --ignore-incomplete
```

## 说明

1. 艾吉泰康UMI
`IGT_UMI_Adapter_and_UDI_Primer` 是 64*64=4096 种 UMI 组合, UMI 长度为 6bp
