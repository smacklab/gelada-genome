#!/bin/sh

mkdir -p id_lists/orthogroups

cp orthogroups/Results_*/Orthogroups/Orthogroups.tsv id_lists/orthogroups/orthogroups.tsv

cat orthogroups/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | sed "s/$/\n/g" | grep -v '^$' > id_lists/orthogroups/sco.txt

mkdir -p palo-tmp

scripts/prep_palo_orthogroups.R

mkdir -p fasta/all_species 
cat fasta/aa/*.faa > fasta/all_species/proteins.faa
cat fasta/cds/*.fna > fasta/all_species/cds.fna

module load contrib/samtools/1.9

samtools faidx fasta/all_species/proteins.faa
samtools faidx fasta/all_species/cds.fna

exit
