#!/bin/bash

source config/settings.cfg

# Use hacked maffilter
export PATH=${HOME}/${bpp_name}/bin:$PATH
export LD_LIBRARY_PATH=${HOME}/${bpp_name}/lib64:$LD_LIBRARY_PATH

# Load mafTools
export PATH=${HOME}/${maftools_name}/bin:$PATH

window_size=$1
perc_intact=$2
perc_identity=$3
n_cores=$4

mkdir -p maf/mrg

# Concatenate all files for the given settings
cat \
	<(echo '##maf version=1') \
	<(cat maf/str/*.restrand.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf | grep -v '^#') > \
	maf/mrg/${maf_version}.restrand.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf

grep "^s $test_species" maf/mrg/${maf_version}.restrand.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf | \
	cut -d ' ' -f 2 | sort -u | sed 's/'${test_species}'\.//g' > \
	data/${test_species}-chr.txt

mkdir -p maf/chr
# Extract chromosomes
cat data/${test_species}-chr.txt | parallel -j $n_cores 'mafExtractor --maf maf/mrg/'${maf_version}'.restrand.window'${window_size}'.intact'${perc_intact}'.identity'${perc_identity}'.maf --seq '${test_species}'.{} --start 0 --stop 9999999999 --soft > maf/chr/'${maf_version}'.{}.window'${window_size}'.intact'${perc_intact}'.identity'${perc_identity}'.maf'

mkdir -p maf/fin
# Sort each chromosome by coordinates
cat data/${test_species}-chr.txt | parallel -j $n_cores 'mafSorter --maf maf/chr/'${maf_version}'.{}.window'${window_size}'.intact'${perc_intact}'.identity'${perc_identity}'.maf --seq '${test_species}'.{} > maf/fin/'${maf_version}'.{}.window'${window_size}'.intact'${perc_intact}'.identity'${perc_identity}'.maf'

mkdir -p maf/com
# Combine into one super maf
cat <(echo '##maf version=1') <(echo '#') \
	<(cat maf/fin/${maf_version}.*.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf | grep -v '^#') > \
	maf/com/${maf_version}.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf

exit