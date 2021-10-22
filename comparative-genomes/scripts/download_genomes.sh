#!/bin/sh

wd=$(pwd)

v=97 # ENSEMBL version

mkdir -p fasta
mkdir -p features

if [ ! -f "fasta/.gitignore" ]
then
	echo '# Ignore everything in this directory'$'\n''*'$'\n''# Except this file'$'\n''!.gitignore' > fasta/.gitignore
fi

if [ ! -f "features/.gitignore" ]
then
	echo '# Ignore everything in this directory'$'\n''*'$'\n''# Except this file'$'\n''!.gitignore' > features/.gitignore
fi

# mkdir -p fasta/aa
# mkdir -p fasta/na
# mkdir -p fasta/cds
# mkdir -p features/gff

for i in $(seq $(wc -l assembly_list.txt | cut -d ' ' -f 1)); do
	acc=$(sed -n ${i}p assembly_list.txt | cut -f 1 | sed "s/ /_/g")
	nam=$(sed -n ${i}p assembly_list.txt | cut -f 2 | sed "s/ /_/g")
	db=$(sed -n ${i}p assembly_list.txt | cut -f 3)
	
	if [ $db = 'refseq' ]; then
		base_url=$(echo ftp://ftp.ncbi.nlm.nih.gov/genomes/all/$(echo $acc | sed "s/\([A-Z]\{3\}\)_\([0-9]\{3\}\)\([0-9]\{3\}\)\([0-9]\{3\}\).*/\1\/\2\/\3\/\4/g")/${acc}_${nam}/${acc}_${nam})
		wget -O fasta/aa/${acc}_${nam}_protein.faa.gz ${base_url}_protein.faa.gz
		wget -O fasta/na/${acc}_${nam}_genomic.fna.gz ${base_url}_genomic.fna.gz
		wget -O fasta/cds/${acc}_${nam}_cds.fna.gz ${base_url}_cds_from_genomic.fna.gz
		wget -O features/gff/${acc}_${nam}_genomic.gff.gz ${base_url}_genomic.gff.gz
	elif [ $db = 'ensembl' ]; then
		base_url=$(echo ftp://ftp.ensembl.org/pub/release-${v})
		wget -O fasta/aa/${acc}_${nam}_protein.faa.gz ${base_url}/fasta/${acc,,}/pep/${acc}.${nam}.pep.all.fa.gz
		wget -O fasta/na/${acc}_${nam}_genomic.fna.gz ${base_url}/fasta/${acc,,}/dna/${acc}.${nam}.dna.toplevel.fa.gz
		wget -O fasta/cds/${acc}_${nam}_cds.fna.gz ${base_url}/fasta/${acc,,}/cds/${acc}.${nam}.cds.all.fa.gz
		wget -O features/gff/${acc}_${nam}_genomic.gff.gz ${base_url}/gff3/${acc,,}/${acc}.${nam}.${v}.gff3.gz
	fi
done;

for i in `ls */*/*.gz`; do gunzip $i; done

module load contrib/samtools/1.9

cd ${wd}/fasta/aa
for i in `ls *.faa`; do if [ ! -f "$i.fai" ]; then samtools faidx $i; fi; done

cd ${wd}/fasta/na
for i in `ls *.fna`; do if [ ! -f "$i.fai" ]; then samtools faidx $i; fi; done

cd ${wd}/fasta/cds
for i in `ls *.fna`; do if [ ! -f "$i.fai" ]; then samtools faidx $i; fi; done

exit
