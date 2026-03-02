# 20s
pigz -dc SURFS15126-RD-SURF5000-T3113V1-250904-HH05-20251125-ZQJ03_R1.fastq.gz > raw.r1.fq &
pigz -dc SURFS15126-RD-SURF5000-T3113V1-250904-HH05-20251125-ZQJ03_R2.fastq.gz > raw.r2.fq &

python fastq_umi_extract-script1.py raw.r1.fq raw.r2.fq

# 20s
pigz filtered.r1.fq
pigz filtered.r2.fq
