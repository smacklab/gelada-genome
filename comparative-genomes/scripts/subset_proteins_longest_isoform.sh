#!/bin/sh

i=$1

mkdir -p id_lists
mkdir -p id_lists/feature_joins
mkdir -p id_lists/aa_longest_isoform

acc=$(sed -n ${i}p assembly_list.txt | cut -f 1 | sed "s/ /_/g")
nam=$(sed -n ${i}p assembly_list.txt | cut -f 2 | sed "s/ /_/g")
db=$(sed -n ${i}p assembly_list.txt | cut -f 3)

if [ ! $db = ensembl ]; then
	cut -f 9 features/gff/${acc}_${nam}_genomic.gff | \
		grep protein_id | grep GeneID | \
		sed "s/GeneID:\([0-9]*\)/\t\1\t/g" | \
		sed "s/Name=\([0-9A-Z_.]*\)/\t\1\t/g" | \
		cut -f 2,4 | sort -u > \
	id_lists/feature_joins/${acc}_${nam}_prot2gene.txt

	grep '>' fasta/cds/${acc}_${nam}_cds.fna | \
		sed 's/>\([A-Za-z0-9|_.]*\) /\1\t/g' | \
		sed 's/protein_id=\([A-Z0-9_.]*\)/\t\1\t/g' | \
		cut -f 1,3 | grep $'\t' | sort -u | \
		awk '{ print $2 "\t" $1 }' > \
	id_lists/feature_joins/${acc}_${nam}_prot2cds.txt
fi

module load contrib/r/3.6.1

scripts/identify_longest_isoform.R ${acc} ${nam} ${db}

module load contrib/samtools/1.9

mkdir -p fasta/aa_longest_isoform

samtools faidx fasta/aa/${acc}_${nam}_protein.faa \
	-r id_lists/aa_longest_isoform/${acc}_${nam}_protein.faa.txt > \
	fasta/aa_longest_isoform/${acc}_${nam}_protein.faa

exit
