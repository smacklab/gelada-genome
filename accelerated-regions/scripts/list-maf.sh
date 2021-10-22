#!/bin/bash

source config/settings.cfg

maf_ftp=${ensembl_root}/maf/ensembl-compara/multiple_alignments/${maf_version}/

mkdir -p data

curl -l $maf_ftp | grep '\.maf\.gz$' > data/maf-files.txt

# Save checksum file
wget -O data/${maf_version}_md5sum.txt ${maf_ftp}MD5SUM

exit