#!/bin/bash

source scripts/_include_options.sh

i=$1
int=$(printf %04d $i)
region=$(sed -n ${i}p data/${genome}_regions.txt)
chr=$(printf %02d $(grep -n $(echo $region | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))
s=$2

if [ ! -f vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.flt.step${s}.vcf.gz.tbi ]; then

tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)

gatk --java-options "-Xmx4g" VariantFiltration \
	--reference ${genome_path} \
	--output  ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.flt.step${s}.vcf.gz \
	--variant  vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz \
	--filter-name "QD" \
	--filter "QD < 2.0" \
	--filter-name "MQ" \
	--filter "MQ < 40.0" \
	--filter-name "FS" \
	--filter "FS > 60.0" \
	--filter-name "MQRS" \
	--filter "MQRankSum < -12.5" \
	--filter-name "RPRS" \
	--filter "ReadPosRankSum < -8.0" \
	--filter-name "SOR" \
	--filter "SOR > 3.0" \
	--missing-values-evaluate-as-failing

if [ $(bgzip -t ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.flt.step${s}.vcf.gz && echo 0 || echo 1) -eq 0 ]; then
mv ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.flt.step${s}.vcf.gz vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.flt.step${s}.vcf.gz
mv ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.flt.step${s}.vcf.gz.tbi vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.flt.step${s}.vcf.gz.tbi
fi

rm -rf $tmp_dir

fi

if [ ! -f vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz.tbi ]; then

tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)

bcftools view -f .,PASS -Oz -o ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.flt.step${s}.vcf.gz

tabix -p vcf ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz

if [ $(bgzip -t ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz && echo 0 || echo 1) -eq 0 ]; then
mv ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz
mv ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz.tbi vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz.tbi
fi

rm -rf $tmp_dir

fi
