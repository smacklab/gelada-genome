#!/bin/bash

i=$1
s=$2

source scripts/_include_options.sh

mkdir -p vcf_chr

bcftools concat -Oz -o vcf_chr/${dataset}.${genome}.bootstrap.chr${i}.step${s}.vcf.gz vcf/${dataset}.${genome}.bootstrap.region.*.chr${i}.flt.step${s}.vcf.gz
