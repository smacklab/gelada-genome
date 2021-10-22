#!/bin/bash

source config/settings.cfg

# Use hacked maffilter
export PATH=${HOME}/${bpp_name}/bin:$PATH
export LD_LIBRARY_PATH=${HOME}/${bpp_name}/lib64:$LD_LIBRARY_PATH

maf=$(sed -n ${1}p data/maf-files.txt)
prefix=$(echo $maf | sed 's/\.maf\.gz$//g')

anchor_species_missing=0

# Add a check to make sure the maf file contains all species
for i in $(echo $anchor_species | sed 's/,/ /g'); do
if zcat maf/raw/${maf} | grep -q '^s '${i}'\.'; then
anchor_species_missing=$(( $anchor_species_missing + 0 ))
else
anchor_species_missing=$(( $anchor_species_missing + 1 ))
fi
done

if zcat maf/raw/${maf} | grep -q '^s '${test_species}'\.'; then
test_species_missing=0
else
test_species_missing=1
fi

if [ $(( $anchor_species_missing + $test_species_missing )) -gt 0 ]; then

echo 'At least one critical species missing from alignment. Skipping.'

echo $maf >> data/maf-fail.txt

else

echo $maf >> data/maf-pass.txt

fi

exit