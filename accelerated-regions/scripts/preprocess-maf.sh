#!/bin/bash

source config/settings.cfg

# Use hacked maffilter
export PATH=${HOME}/${bpp_name}/bin:$PATH
export LD_LIBRARY_PATH=${HOME}/${bpp_name}/lib64:$LD_LIBRARY_PATH

# Load mafTools
export PATH=${HOME}/${maftools_name}/bin:$PATH

maf=$(sed -n ${1}p data/maf-pass.txt)
prefix=$(echo $maf | sed 's/\.maf\.gz$//g')
species_order=$(echo $all_species | sed 's/^\(.*\),'${ref_species}'\(.*\)/'${ref_species}',\1\2/g')

mkdir -p maf/pre

# Subset species
maffilter input.file=maf/raw/${maf} input.file.compression=gzip input.dots=error maf.filter='
Subset(
	species=('${all_species}'),
	strict=no,
	remove_duplicates=no),
Output(
	file='maf/pre/${prefix}.species.maf',
	compression=none,mask=yes)'

# Remove duplicate species
mafDuplicateFilter --maf maf/pre/${prefix}.species.maf > \
	maf/pre/${prefix}.rmdup.maf

# Move reference species to first row of block
mafRowOrderer \
	--maf maf/pre/${prefix}.rmdup.maf \
	--order ${species_order} > \
	maf/pre/${prefix}.reference.maf

# Reference all blocks to + strand on ref_species
mafStrander --maf maf/pre/${prefix}.reference.maf \
	--seq ${ref_species} --strand + > \
	maf/pre/${prefix}.stranded.maf

# Sort blocks by genomic position on reference
mafSorter --maf maf/pre/${prefix}.stranded.maf \
	--seq ${ref_species} > \
	maf/pre/${prefix}.sorted.maf

# Also create gzipped version
cat maf/pre/${prefix}.sorted.maf | gzip > \
	maf/pre/${prefix}.sorted.maf.gz

exit