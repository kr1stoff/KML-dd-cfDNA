/home/mengxf/miniforge3/envs/basic2/bin/picard FastqToSam \
	F1=/data/mengxf/Task/KML260115-dd-cfDNA-iGT/results/260126/qc/fastp/ZQJ12.1.fastq.gz \
	F2=/data/mengxf/Task/KML260115-dd-cfDNA-iGT/results/260126/qc/fastp/ZQJ12.2.fastq.gz \
	O=ZQJ12.unmapped.bam \
	SAMPLE_NAME=ZQJ12

/home/mengxf/miniforge3/envs/basic2/bin/fgbio ExtractUmisFromBam \
	--input=ZQJ12.unmapped.bam --output=ZQJ12.unmapped.withUMI.bam \
	--read-structure=6M144T 6M144T --molecular-index-tags=ZA ZB \
	--single-tag=RX

/home/mengxf/miniforge3/envs/basic2/bin/picard SamToFastq I=ZQJ12.unmapped.withUMI.bam F=ZQJ12.SamToFastq INTERLEAVE=true

bwa mem -p -t 32 /data/mengxf/Database/reference/hg19/hg19.fa ZQJ12.SamToFastq >ZQJ12.mapped.bam

/home/mengxf/miniforge3/envs/basic2/bin/picard MergeBamAlignment \
	UNMAPPED=ZQJ12.unmapped.withUMI.bam \
	ALIGNED=ZQJ12.mapped.bam \
	O=ZQJ12.mapped.withUMI.bam \
	R=/data/mengxf/Database/reference/hg19/hg19.fa \
	SO=coordinate \
	ALIGNER_PROPER_PAIR_FLAGS=true \
	MAX_GAPS=-1 \
	ORIENTATIONS=FR \
	VALIDATION_STRINGENCY=SILENT \
	CREATE_INDEX=true

/home/mengxf/miniforge3/envs/basic2/bin/fgbio GroupReadsByUmi \
	--input=ZQJ12.mapped.withUMI.bam --output=ZQJ12.GroupedReads.bam \
	--strategy=paired --edits=1 --min-map-q=20

/home/mengxf/miniforge3/envs/basic2/bin/fgbio CallDuplexConsensusReads \
	--input=ZQJ12.GroupedReads.bam --output=ZQJ12.consensus.unmapped.bam \
	--error-rate-pre-umi=45 --error-rate-post-umi=30 \
	--min-input-base-quality=30

/home/mengxf/miniforge3/envs/basic2/bin/fgbio FilterConsensusReads -Xmx8G \
	--input=ZQJ12.consensus.unmapped.bam --output=ZQJ12.consensus.filtered.unmapped.bam \
	--ref=/data/mengxf/Database/reference/hg19/hg19.fa \
	--min-reads=2 1 1 \
	--max-read-error-rate=0.05 \
	--max-base-error-rate=0.1 \
	--min-base-quality=50 \
	--max-no-calls=0.05

/home/mengxf/miniforge3/envs/basic2/bin/picard SamToFastq I=ZQJ12.consensus.filtered.unmapped.bam F=ZQJ12.consensus.filtered.unmapped.Fastq INTERLEAVE=true

bwa mem -p -t 32 /data/mengxf/Database/reference/hg19/hg19.fa ZQJ12.consensus.filtered.unmapped.Fastq >ZQJ12.consensus.filtered.mapped.bam

/home/mengxf/miniforge3/envs/basic2/bin/picard SortSam INPUT=ZQJ12.consensus.filtered.mapped.bam OUTPUT=ZQJ12.consensus.filtered.mapped.sorted.bam SORT_ORDER=queryname

/home/mengxf/miniforge3/envs/basic2/bin/picard SortSam -Xmx8G INPUT=ZQJ12.consensus.filtered.unmapped.bam OUTPUT=ZQJ12.consensus.filtered.unmapped.sorted.bam SORT_ORDER=queryname

/home/mengxf/miniforge3/envs/basic2/bin/picard MergeBamAlignment -Xmx8G \
	UNMAPPED=ZQJ12.consensus.filtered.unmapped.sorted.bam \
	ALIGNED=ZQJ12.consensus.filtered.mapped.sorted.bam \
	O=ZQJ12.consensus.filtered.mapped.withUMI.bam \
	R=/data/mengxf/Database/reference/hg19/hg19.fa \
	SO=coordinate \
	ALIGNER_PROPER_PAIR_FLAGS=true \
	MAX_GAPS=-1 \
	ORIENTATIONS=FR \
	VALIDATION_STRINGENCY=SILENT \
	CREATE_INDEX=true

samtools sort -n -u ZQJ12.consensus.filtered.mapped.withUMI.bam |
	/home/mengxf/miniforge3/envs/basic2/bin/fgbio ClipBam -Xmx8G \
		--input=/dev/stdin \
		--output=ZQJ12.consensus.filtered.mapped.withUMI.clipped.bam \
		--ref=/data/mengxf/Database/reference/hg19/hg19.fa \
		--clipping-mode=Soft \
		--clip-overlapping-reads=true

# * part 20260130
# 排序, 索引
samtools sort --threads 8 \
	ZQJ12.consensus.filtered.mapped.withUMI.clipped.bam \
	-o ZQJ12.consensus.filtered.mapped.withUMI.clipped.sorted.bam

samtools index ZQJ12.consensus.filtered.mapped.withUMI.clipped.sorted.bam

# mpileup → call → bgzip, 索引
bcftools mpileup \
	--threads 8 \
	--max-depth 25000 \
	--min-MQ 20 \
	--min-BQ 30 \
	--no-BAQ \
	-Ou \
	-f /data/mengxf/Database/reference/hg19/hg19.fa \
	-R /data/mengxf/GitHub/KML-dd-cfDNA/wf-igt-500snp/assets/probeCov.predictSlop500bp.bed \
	ZQJ12.sorted.bam |
	bcftools call \
		--threads 8 \
		--multiallelic-caller \
		-Ov |
	bgzip -c \
		>ZQJ12.region.call.vcf.gz

bcftools index ZQJ12.region.call.vcf.gz

# 提取 500 SNP
bcftools view \
	-R /data/mengxf/Task/KML260115-dd-cfDNA-iGT/Chimerism_SNP_Panel-T3113V1/T3113V1/anno/loci.anno.vcf.gz \
	ZQJ12.region.call.vcf.gz |
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%AD]\n' >ZQJ12.500snps.txt

# 计算 500 SNP 的平均深度
awk '$5 != "." {sum+=$5; n++} END {if(n>0) printf int(sum/n)}' ZQJ12.500snps.txt >ZQJ12.500snps.avgdepth.txt

# 提取 2000 DP 以上, 不在目标区域的 SNP
bcftools view \
	-i 'ALT!="." && DP>=2000 && DP!="."' \
	ZQJ12.region.call.vcf.gz |
	bedtools intersect -header -v -wa \
		-a - \
		-b /data/mengxf/Task/KML260115-dd-cfDNA-iGT/Chimerism_SNP_Panel-T3113V1/T3113V1/loci.bed |
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
		>ZQJ12.region.2000dp.no-target.txt

# 计算测序错误率
# 计算 dd-cfDNA 比率
