#!/bin/bash

source config/settings.cfg

# Use hacked maffilter
export PATH=${HOME}/${bpp_name}/bin:$PATH
export LD_LIBRARY_PATH=${HOME}/${bpp_name}/lib64:$LD_LIBRARY_PATH

maf=$(sed -n ${1}p data/maf-pass.txt)
prefix=$(echo $maf | sed 's/\.maf\.gz$//g')

window_size=$2
perc_intact=$3

# Maximum number of gapped positions
max_gap=$(echo "$window_size-$window_size*$perc_intact/100" | bc -l | awk '{printf "%d\n",$0+=$0<0?-1:0}')

mkdir -p maf/aln

mkdir -p mafft/scripts
mkdir -p mafft/fasta
mkdir -p mafft/logs

cat data/run-mafft.sh | \
	sed 's:blockIn:mafft/fasta/'${prefix}'_in:g' | \
	sed 's:blockOut:mafft/fasta/'${prefix}'_out:g' | \
	sed 's:mafft.log:mafft/logs/'${prefix}'_log.log:g'> \
	mafft/scripts/run_mafft_${prefix}.sh

chmod +x mafft/scripts/run_mafft_${prefix}.sh

# Remove gaps
maffilter input.file=maf/pre/${prefix}.sorted.maf.gz input.file.compression=gzip input.dots=error maf.filter='
Subset(
	species=('${anchor_species}','${test_species}'),
	strict=yes,
	keep=no,
	remove_duplicates=yes),
XFullGap(
	species=('${anchor_species}','${test_species}')),
SystemCall(
	name=MAFFT,
	input.file=mafft/fasta/'${prefix}'_in.fasta,
	input.format=Fasta,
	output.file=mafft/fasta/'${prefix}'_out.fasta,
	output.format=Fasta,
	call=mafft/scripts/run_mafft_'${prefix}'.sh),
MinBlockLength(
	min_length='${window_size}'),
AlnFilter(
	species=('${anchor_species}'),
	window.size='${window_size}',
	window.step=1,
	max.ent=-1,
	max.gap='${max_gap}',
	missing_as_gap=yes,
	file='maf/aln/${prefix}.alnfilter.window${window_size}.intact${perc_intact}.maf.gz',
	compression=gzip)'

exit