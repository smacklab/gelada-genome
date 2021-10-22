#!/bin/bash

i=$1
j=$2

source scripts/_include_options.sh

int=$(printf %04d $(grep -n $j data/${genome}_regions.txt | cut -d ':' -f 1))
chr=$(printf %02d $(grep -n $(echo $j | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))

tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)

java -jar $PICARDJAR AddOrReplaceReadGroups \
	I=bam/${i}.${genome}.region.${int}.chr${chr}.bam \
	O=${tmp_dir}/${i}.${genome}.region.${int}.chr${chr}.bam \
	CREATE_INDEX=true \
	RGID=$i RGLB=$i RGPL=Illumina RGSM=$i RGPU=1

if [ $(samtools quickcheck ${tmp_dir}/${i}.${genome}.region.${int}.chr${chr}.bam && echo 0 || echo 1) -eq 0 ]; then

mv ${tmp_dir}/${i}.${genome}.region.${int}.chr${chr}.bam bam/${i}.${genome}.region.${int}.chr${chr}.bam
mv ${tmp_dir}/${i}.${genome}.region.${int}.chr${chr}.bai bam/${i}.${genome}.region.${int}.chr${chr}.bam.bai

rm -rf $tmp_dir

fi

