#!/bin/bash

source scripts/_include_options.sh

i=$1
int=$(printf %04d $i)
region=$(sed -n ${i}p data/${genome}_regions.txt)
chr=$(printf %02d $(grep -n $(echo $region | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))
s=$2

# vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz

mkdir -p vcf_snps

vcftools \
	--gzvcf vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz \
	--recode \
	--remove-indels \
	--maf 0.01 \
	--max-missing 0.9 \
	--stdout | \
	bgzip -c > vcf_snps/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.good.step${s}.vcf.gz

tabix -p vcf vcf_snps/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.good.step${s}.vcf.gz

# vcf-concat vcf_snps/gelada.tgel1.bootstrap.region.*.chr{01..21}.good.step2.vcf.gz | bgzip -c > vcf_final/gelada.tgel1.filtered.all.step2.vcf.gz
# vcf-concat vcf_snps/gelada.panubis1.bootstrap.region.*.chr{01..20}.good.step2.vcf.gz | bgzip -c > vcf_final/gelada.panubis1.filtered.all.step2.vcf.gz
