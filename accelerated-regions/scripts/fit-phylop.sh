#!/bin/bash

source config/settings.cfg

# Load phast
export PATH=${HOME}/${phast_name}/bin:$PATH

window_size=$2
perc_intact=$3
perc_identity=$4

mkdir -p mod/ind

maf=$(sed -n ${1}p data/${test_species}-blocks.txt)

# Extract chromosome
chr=$(echo $maf | sed 's:maf/ind/'${maf_version}'\.\(.*\)\.window'${window_size}'\.intact'${perc_intact}'\.identity'${perc_identity}'\.\([0-9]*\)\.maf:\1:g')
# Extract block number
block=$(echo $maf | sed 's:maf/ind/'${maf_version}'\.\(.*\)\.window'${window_size}'\.intact'${perc_intact}'\.identity'${perc_identity}'\.\([0-9]*\)\.maf:\2:g')
# Extract block starting position
pos=$(grep 's '${test_species}.${chr} $maf | tr -s ' ' | cut -d ' ' -f 3)
# Extract block length
len=$(grep 's '${test_species}.${chr} $maf | tr -s ' ' | cut -d ' ' -f 4)
# Extract strand (should be + always)
str=$(grep 's '${test_species}.${chr} $maf | tr -s ' ' | cut -d ' ' -f 5)

mkdir -p mod/acc

phyloFit \
	--tree "$tree" \
	--subst-mod REV \
	--msa-format MAF \
	--out-root mod/acc/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block} \
	$maf

mkdir -p results/phylop

phyloP \
	--method LRT \
	--mode CONACC \
	--msa-format MAF \
	--base-by-base \
	--branch $test_species \
	mod/acc/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block}.mod \
	$maf > \
	results/phylop/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block}.wig

# mkdir -p results/acc
# 
# csplit -z -s -b '%04d.txt' \
# 	-f results/acc/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}. \
# 	results/phylop/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block}.wig \
# 	'/^fixedStep/' '{*}'
# 
# # Delete first file, which is just the header
# #null_scale alt_scale alt_subscale lnlratio pval
# 
# rm -rf results/acc/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.0000.txt

exit