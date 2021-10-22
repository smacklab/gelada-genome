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

phyloFit \
	--subst-mod REV \
	--init-model mod/all/${maf_version}.window${window_size}.intact${perc_intact}.identity${perc_identity}.mod \
	--scale-only \
	--msa-format MAF \
	--out-root mod/ind/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block}.nul \
	$maf

phyloFit \
	--subst-mod REV \
	--init-model mod/all/${maf_version}.window${window_size}.intact${perc_intact}.identity${perc_identity}.mod \
	--scale-only \
	--scale-subtree ${test_species}:loss \
	--msa-format MAF \
	--out-root mod/ind/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block}.alt \
	$maf

lnl_nul=$(grep '^TRAINING_LNL: ' mod/ind/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block}.nul.mod | cut -d ':' -f 2 | xargs)
lnl_alt=$(grep '^TRAINING_LNL: ' mod/ind/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block}.alt.mod | cut -d ':' -f 2 | xargs)

# Save results
mkdir -p results/cons

echo ${chr}_${block}$'\t'$chr$'\t'$pos$'\t'$str$'\t'$len$'\t'$lnl_nul$'\t'$lnl_alt > \
	results/cons/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.${block}.lrt.txt

exit