#!/bin/bash

source config/settings.cfg

chr=$(sed -n ${1}p data/${test_species}-chr.txt)

window_size=$2
perc_intact=$3
perc_identity=$4

mkdir -p maf/ind

csplit -z -s -b '%06d.maf' \
	-f maf/ind/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}. \
	maf/fin/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.maf \
	'/^a/' '{*}'

rm -rf maf/ind/${maf_version}.${chr}.window${window_size}.intact${perc_intact}.identity${perc_identity}.000000.maf

exit