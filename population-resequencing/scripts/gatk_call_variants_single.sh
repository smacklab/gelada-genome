#!/bin/bash

i=$1
j=$2
s=$3

source scripts/_include_options.sh

int=$(printf %04d $(grep -n $j data/${genome}_regions.txt | cut -d ':' -f 1))
chr=$(printf %02d $(grep -n $(echo $j | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))

if [ $s -gt 0 ]; then
suffix=.step${s}
else
suffix=
fi

module load htslib/1.9.0

mkdir -p gvcf

if [ ! -f gvcf/${i}.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz.tbi ]; then

tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)

gatk --java-options '-Xmx4g' HaplotypeCaller \
	--reference ${genome_path} \
	--output ${tmp_dir}/${i}.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz \
	--input bam/${i}.${genome}.region.${int}.chr${chr}${suffix}.bam \
	--intervals $j --emit-ref-confidence GVCF

if [ $(bgzip -t ${tmp_dir}/${i}.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz && echo 0 || echo 1) -eq 0 ]; then
mv ${tmp_dir}/${i}.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz gvcf/${i}.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz
mv ${tmp_dir}/${i}.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz.tbi gvcf/${i}.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz.tbi
fi

rm -rf $tmp_dir

fi
