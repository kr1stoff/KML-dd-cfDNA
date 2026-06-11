POSITIONS="/data/mengxf/Task/KML260115-dd-cfDNA-iGT/work/260512-compare-process/positions.txt"
LOCI_VCF="/data/mengxf/Task/KML260115-dd-cfDNA-iGT/work/260512-compare-process/loci.anno.vcf"
# 3. raw + clipOverlap
/home/mengxf/miniforge3/envs/basic2/bin/bam clipOverlap --in ${BAM} --out ${SAMPLE}.clipOverlap.bam --stats
samtools index ${SAMPLE}.clipOverlap.bam
cat ${POSITIONS} | while read position;do echo "/home/mengxf/miniforge3/envs/basic/bin/samtools mpileup --fasta-ref ${REFERENCE} --max-depth 500000 --min-BQ 30 -r ${position} ${SAMPLE}.clipOverlap.bam > tmp/${position}.mpileup";done | parallel -j 8
cat tmp/*.mpileup > ${SAMPLE}.clipOverlap.mpileup
rm tmp/*.mpileup
python /data/mengxf/Task/KML260115-dd-cfDNA-iGT/work/260512-compare-process/parse_bases.py ${SAMPLE}.clipOverlap.mpileup ${SAMPLE}.clipOverlap.base_stat.tsv
python /data/mengxf/Task/KML260115-dd-cfDNA-iGT/work/260512-compare-process/calc_af.py ${LOCI_VCF} ${SAMPLE}.clipOverlap.base_stat.tsv ${SAMPLE}.clipOverlap.alt_freq.tsv
