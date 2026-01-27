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

bwa mem -p -t 32 /data/mengxf/Database/reference/hg19/hg19.fa ZQJ12.SamToFastq > ZQJ12.mapped.bam

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

# * 限制内存
/home/mengxf/miniforge3/envs/basic2/bin/fgbio FilterConsensusReads -Xmx8G \
    --input=ZQJ12.consensus.unmapped.bam --output=ZQJ12.consensus.filtered.unmapped.bam \
    --ref=/data/mengxf/Database/reference/hg19/hg19.fa \
    --min-reads=2 1 1 \
    --max-read-error-rate=0.05 \
    --max-base-error-rate=0.1 \
    --min-base-quality=50 \
    --max-no-calls=0.05

/home/mengxf/miniforge3/envs/basic2/bin/picard SamToFastq I=ZQJ12.consensus.filtered.unmapped.bam F=ZQJ12.consensus.filtered.unmapped.Fastq INTERLEAVE=true

bwa mem -p -t 32 /data/mengxf/Database/reference/hg19/hg19.fa ZQJ12.consensus.filtered.unmapped.Fastq > ZQJ12.consensus.filtered.mapped.bam

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

# /home/mengxf/miniforge3/envs/basic2/bin/fgbio ClipBam -Xmx8G \
#     --input=ZQJ12.consensus.filtered.mapped.withUMI.bam --output=ZQJ12.consensus.filtered.mapped.withUMI.clipped.bam \
#     --ref=/data/mengxf/Database/reference/hg19/hg19.fa --clipping-mode=Soft --clip-overlapping-reads=true

/home/mengxf/miniforge3/envs/basic2/bin/samtools sort -n -u ZQJ12.consensus.filtered.mapped.withUMI.bam | \
    /home/mengxf/miniforge3/envs/basic2/bin/fgbio ClipBam -Xmx8G \
    --input=/dev/stdin \
    --output=ZQJ12.consensus.filtered.mapped.withUMI.clipped.bam \
    --ref=/data/mengxf/Database/reference/hg19/hg19.fa \
    --clipping-mode=Soft \
    --clip-overlapping-reads=true
