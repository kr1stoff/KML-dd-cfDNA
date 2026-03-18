# 20s
time pigz -dc SURFS15126-RD-SURF5000-T3113V1-250904-HH05-20251125-ZQJ03_R1.fastq.gz > raw.r1.fq &
time pigz -dc SURFS15126-RD-SURF5000-T3113V1-250904-HH05-20251125-ZQJ03_R2.fastq.gz > raw.r2.fq

# 1min
time /data/mengxf/GitHub/KML-dd-cfDNA/rust/umi_filter/target/release/umi_filter raw.r1.fq raw.r2.fq rsout.1.fq rsout.2.fq /data/mengxf/GitHub/KML-dd-cfDNA/wf-igt-500snp/assets/umi.txt rsumi_stats.txt

# 20s
time pigz -c rsout.1.fq > rsout.1.fq.gz &
time pigz -c rsout.2.fq > rsout.2.fq.gz
