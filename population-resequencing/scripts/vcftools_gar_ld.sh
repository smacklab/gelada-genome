#!/bin/bash

mkdir -p gar
cd gar

zcat ../vcf/gelada.tgel1.bootstrap.region.0866.chr20.pas.step2.vcf.gz > ZNF536_gelada.vcf
zcat ../vcf/gelada.tgel1.bootstrap.region.0895.chr21.pas.step2.vcf.gz > RBFOX1_gelada.vcf

vcftools --recode --remove-indels --min-alleles 2 --max-alleles 2 --vcf RBFOX1_gelada.vcf --stdout | bgzip -c > RBFOX1_gelada_biallelic.vcf.gz
vcftools --recode --remove-indels --min-alleles 2 --max-alleles 2 --vcf ZNF536_gelada.vcf --stdout | bgzip -c > ZNF536_gelada_biallelic.vcf.gz

tabix -p vcf RBFOX1_gelada_biallelic.vcf.gz
tabix -p vcf ZNF536_gelada_biallelic.vcf.gz

echo 'NC_037688.1'$'\t''75989659'$'\t''75989760'$'\n''NC_037688.1'$'\t''76699130'$'\t''76699287' > RBFOX1_regions.txt
echo 'NC_037687.1'$'\t''49184834'$'\t''49184933'$'\n''NC_037687.1'$'\t''49695969'$'\t''49696366' > ZNF536_regions.txt

tabix -h -R RBFOX1_regions.txt RBFOX1_gelada_biallelic.vcf.gz > RBFOX1_gelada_biallelic_region.vcf
tabix -h -R ZNF536_regions.txt ZNF536_gelada_biallelic.vcf.gz > ZNF536_gelada_biallelic_region.vcf

vcftools --vcf RBFOX1_gelada_biallelic_region.vcf --geno-r2 --out foo --ld-window-bp-min 500000 --out RBFOX1_ld
vcftools --vcf ZNF536_gelada_biallelic_region.vcf --geno-r2 --out foo --ld-window-bp-min 500000 --out ZNF536_ld

# Add 1000 bp on each end of each GAR
echo 'NC_037688.1'$'\t''75988659'$'\t''75990760'$'\n''NC_037688.1'$'\t''76698130'$'\t''76700287' > RBFOX1_regions2.txt
echo 'NC_037687.1'$'\t''49183834'$'\t''49185933'$'\n''NC_037687.1'$'\t''49694969'$'\t''49697366' > ZNF536_regions2.txt

tabix -h -R RBFOX1_regions2.txt RBFOX1_gelada_biallelic.vcf.gz > RBFOX1_gelada_biallelic_region2.vcf
tabix -h -R ZNF536_regions2.txt ZNF536_gelada_biallelic.vcf.gz > ZNF536_gelada_biallelic_region2.vcf

vcftools --vcf RBFOX1_gelada_biallelic_region2.vcf --geno-r2 --out foo --ld-window-bp-min 500000 --out RBFOX1_ld2
vcftools --vcf ZNF536_gelada_biallelic_region2.vcf --geno-r2 --out foo --ld-window-bp-min 500000 --out ZNF536_ld2
