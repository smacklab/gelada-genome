#!/bin/bash

i=$1
j=$2
s=$3

source scripts/_include_options.sh

int=$(printf %04d $(grep -n $j data/${genome}_regions.txt | cut -d ':' -f 1))
chr=$(printf %02d $(grep -n $(echo $j | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))

mkdir -p reports

if [ ! -f reports/${i}.${genome}.region.${int}.chr${chr}.bqr.step${s}.txt ]; then

gatk --java-options '-Xmx4g' BaseRecalibrator \
	--reference ${genome_path} \
	--input bam/${i}.${genome}.region.${int}.chr${chr}.bam \
	--use-original-qualities \
	--output reports/${i}.${genome}.region.${int}.chr${chr}.bqr.step${s}.txt \
	--known-sites vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step$(($s-1)).vcf.gz

fi

if [ -f reports/${i}.${genome}.region.${int}.chr${chr}.bqr.step${s}.txt ]; then

if [ ! -f bam/${i}.${genome}.region.${int}.chr${chr}.step${s}.bam ] || [ $(samtools quickcheck bam/${i}.${genome}.region.${int}.chr${chr}.step${s}.bam && echo 0 || echo 1) -eq 1 ]; then

gatk --java-options '-Xmx4g' ApplyBQSR \
	--add-output-sam-program-record \
	--reference ${genome_path} \
	--input bam/${i}.${genome}.region.${int}.chr${chr}.bam \
	--use-original-qualities \
	--output bam/${i}.${genome}.region.${int}.chr${chr}.step${s}.bam \
	--bqsr-recal-file reports/${i}.${genome}.region.${int}.chr${chr}.bqr.step${s}.txt

mv bam/${i}.${genome}.region.${int}.chr${chr}.step${s}.bai bam/${i}.${genome}.region.${int}.chr${chr}.step${s}.bam.bai

fi

else

echo 'Recalibration file failed to generate'
exit 113

fi
