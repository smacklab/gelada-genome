#!/bin/bash

source config/settings.cfg

# Use hacked maffilter
export PATH=${HOME}/${bpp_name}/bin:$PATH
export LD_LIBRARY_PATH=${HOME}/${bpp_name}/lib64:$LD_LIBRARY_PATH

maf=$(sed -n ${1}p data/maf-pass.txt)
prefix=$(echo $maf | sed 's/\.maf\.gz$//g')

window_size=$2
perc_intact=$3
perc_identity=$4

# Maximum number of variable positions
max_pos=$(echo "$window_size-$window_size*$perc_identity/100" | bc -l | awk '{printf "%d\n",$0+=$0<0?-1:0}')

mkdir -p maf/ent

# Remove variable regions
maffilter input.file=maf/aln/${prefix}.alnfilter.window${window_size}.intact${perc_intact}.maf.gz input.file.compression=gzip maf.filter='
EntropyFilter(
	species=('${anchor_species}'),
	window.size='${window_size}',
	window.step=1,
	max.ent=0.1,
	max.pos='${max_pos}',
	missing_as_gap=yes,
	ignore_gaps=no,
	file='maf/ent/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf.gz',
	compression=gzip)'

exit