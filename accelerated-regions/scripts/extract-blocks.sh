#!/bin/bash

source config/settings.cfg

# Use hacked maffilter
export PATH=${HOME}/${bpp_name}/bin:$PATH
export LD_LIBRARY_PATH=${HOME}/${bpp_name}/lib64:$LD_LIBRARY_PATH

# Load phast
export PATH=${HOME}/${phast_name}/bin:$PATH

# Load mafTools
export PATH=${HOME}/${maftools_name}/bin:$PATH

maf=$(sed -n ${1}p data/maf-pass.txt)
prefix=$(echo $maf | sed 's/\.maf\.gz$//g')

window_size=$2
perc_intact=$3
perc_identity=$4

mkdir -p reports/ent
mkdir -p reports/crd

# Write statistics and output coordinates
maffilter input.file=maf/ent/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf.gz input.file.compression=gzip maf.filter='
SequenceStatistics(
	statistics=(BlockLength,SequenceLength(species='${test_species}')),
    ref_species='${ref_species}',
    file='reports/ent/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.stats.txt',
    compression=none),
OutputCoordinates(
	file='reports/crd/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.coord.txt',
	compression=none,
	species=('${ref_species}'),
	output_src_size=no)'

mkdir -p reports/gff

# Convert to gff3
tail -n +2 reports/crd/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.coord.txt | \
	awk '{print $1"\t.\tblock\t"$3+1"\t"$4"\t.\t"$2"\t.\tID="$3}' > \
	reports/gff/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.gff3

if [ $(wc -l reports/gff/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.gff3 | xargs | cut -d ' ' -f 1) -gt 0 ]; then

mkdir -p maf/ext

# Extract coordinates in GFF3 file (add header because it will be required by mafTools)
cat \
<(echo '##maf version=1') \
<(maf_parse \
	--features reports/gff/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.gff3 \
	maf/pre/${prefix}.sorted.maf) > \
	maf/ext/${prefix}.extract.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf

mkdir -p maf/rln

mkdir -p mafft/scripts
mkdir -p mafft/fasta
mkdir -p mafft/logs

cat data/run-mafft.sh | \
	sed 's:blockIn:mafft/fasta/'${prefix}'_in:g' | \
	sed 's:blockOut:mafft/fasta/'${prefix}'_out:g' | \
	sed 's:mafft.log:mafft/logs/'${prefix}'_log.log:g'> \
	mafft/scripts/run_mafft_${prefix}.sh

chmod +x mafft/scripts/run_mafft_${prefix}.sh

maffilter input.file=maf/ext/${prefix}.extract.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf input.file.compression=none input.dots=error maf.filter='
SystemCall(
	name=MAFFT,
	input.file=mafft/fasta/'${prefix}'_in.fasta,
	input.format=Fasta,
	output.file=mafft/fasta/'${prefix}'_out.fasta,
	output.format=Fasta,
	call=mafft/scripts/run_mafft_'${prefix}'.sh),
Output(
	file='maf/rln/${prefix}.realign.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf',
	compression=none,mask=yes)'

mkdir -p maf/srt

# Sort rows by all_species (test branch should be first)
mafRowOrderer \
	--maf maf/rln/${prefix}.realign.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf \
	--order ${all_species} > \
	maf/srt/${prefix}.reorder.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf

mkdir -p maf/str

# Fix strandedness for test species (mafStrander fails for empty maf files)
mafStrander --maf maf/srt/${prefix}.reorder.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf \
	--seq ${test_species} --strand + > \
	maf/str/${prefix}.restrand.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf

else

echo 'Warning: no blocks found in 'reports/gff/${prefix}.entfilter.window${window_size}.intact${perc_intact}.identity${perc_identity}.gff3
cat <(echo '##maf version=1') <(echo '#') > \
	maf/str/${prefix}.restrand.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf

fi

exit