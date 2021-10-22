#!/bin/bash

source config/settings.cfg

maf_ftp=${ensembl_root}/maf/ensembl-compara/multiple_alignments/${maf_version}/
cmd_prefix="wget $maf_ftp"

maf=$(sed -n ${1}p data/maf-files.txt)

# Construct wget command
cmd=${cmd_prefix}${maf}

mkdir -p maf/raw
cd maf/raw

$cmd

if [ $(grep ' '$maf'$' data/${maf_version}_md5sum.txt | cut -d ' ' -f 1) = $(md5sum maf/raw/${maf} | cut -d ' ' -f 1) ]; then
i=10
else
i=1
fi

# If MD5 check fails, try again up to 10 total tries

while [ $i -lt 10 ]; do
echo 'MD5 check failed. Attempt #'$(( $i + 1 ))'.'

$cmd
if [ $(grep ' '$maf'$' data/${maf_version}_md5sum.txt | cut -d ' ' -f 1) = $(md5sum maf/raw/${maf} | cut -d ' ' -f 1) ]; then
i=10
else
i=$(( $i + 1 ))
fi
done

exit