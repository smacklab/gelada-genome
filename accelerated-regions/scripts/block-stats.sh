#!/bin/bash

source config/settings.cfg

# Use hacked maffilter
export PATH=${HOME}/${bpp_name}/bin:$PATH
export LD_LIBRARY_PATH=${HOME}/${bpp_name}/lib64:$LD_LIBRARY_PATH

maf=$(sed -n ${1}p data/${test_species}-blocks.txt)
prefix=$(echo $maf | sed 's/maf\/ind\/\(.*\)\.maf$/\1/g')

mkdir -p reports/fin

# Write statistics and output coordinates
maffilter input.file=${maf} input.file.compression=none maf.filter='
SequenceStatistics(
	statistics=(BlockSize,BlockLength,BlockCounts,'$(echo 'SequenceLength(species='$(echo $all_species | sed 's/,/),SequenceLength(species=/g')')')'),
    ref_species='${test_species}',
    file='reports/fin/${prefix}.stats.txt',
    compression=none)'

mkdir -p reports/ref

# Write statistics and output coordinates
maffilter input.file=${maf} input.file.compression=none maf.filter='
SequenceStatistics(
	statistics=(BlockSize,BlockLength,BlockCounts,'$(echo 'SequenceLength(species='$(echo $all_species | sed 's/,/),SequenceLength(species=/g')')')'),
    ref_species='${ref_species}',
    file='reports/ref/${prefix}.stats.${ref_species}.txt',
    compression=none)'

exit