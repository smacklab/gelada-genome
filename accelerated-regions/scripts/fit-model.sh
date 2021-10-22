#!/bin/bash

source config/settings.cfg

# Load phast
export PATH=${HOME}/${phast_name}/bin:$PATH

window_size=$1
perc_intact=$2
perc_identity=$3

mkdir -p mod/all

phyloFit \
	--tree "$tree" \
	--subst-mod REV \
	--msa-format MAF \
	--out-root mod/all/${maf_version}.window${window_size}.intact${perc_intact}.identity${perc_identity} \
	maf/com/${maf_version}.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf

exit