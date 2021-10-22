#!/usr/bin/sh

export NCBI_API_KEY=$(cat config/ncbi_api_key.txt)

mkdir -p tmp
mkdir -p tmp/cds

# module load contrib/edirect/12.0
module load edirect/1.0

# for i in $(seq $(wc -l id_lists/orthogroups/sco_all.txt | cut -d ' ' -f 1)); do
for i in $(seq $(wc -l id_lists/orthogroups/sco_most.txt | cut -d ' ' -f 1)); do
	echo $i;
# 	orthogroup=`sed -n ${i}p id_lists/orthogroups/sco_all.txt`
	orthogroup=`sed -n ${i}p id_lists/orthogroups/sco_most.txt`;
	esearch -db protein -query $(cat id_lists/aa_palo_isoform/${orthogroup}.txt| grep -v ^ENS | sed ':a;N;$!ba;s/\n/,/g') | efetch -format fasta_cds_na > tmp/cds/${orthogroup}.fa;
done

cat id_lists/cds_palo_isoform/*.txt | grep '^ENS' > tmp/cds_ensembl.txt

# module load contrib/samtools/1.9
module load samtools/1.9

samtools faidx fasta/all_species/cds.fna -r tmp/cds_ensembl.txt > tmp/cds/ENSEMBL.fa

cat tmp/cds/*.fa > fasta/all_species/cds_corrected.fna

samtools faidx fasta/all_species/cds_corrected.fna

# module load contrib/r/3.6.1
module load r/latest

scripts/fix_cds_accession_lists.R

cp fasta/all_species/cds_corrected.fna fasta/all_species/cds_corrected_all_codons.fna
scripts/fix_cds_codons.R

samtools faidx fasta/all_species/cds_corrected_all_codons.fna

# module load contrib/seqkit/0.10.0
module load seqkit/0.12.0

mkdir -p fasta/sco/cds_corrected

# for i in $(seq $(wc -l id_lists/orthogroups/sco_all.txt | cut -d ' ' -f 1)); do
for i in $(seq $(wc -l id_lists/orthogroups/sco_most.txt | cut -d ' ' -f 1)); do
echo $i;
# orthogroup=`sed -n ${i}p id_lists/orthogroups/sco_all.txt`;
orthogroup=`sed -n ${i}p id_lists/orthogroups/sco_most.txt`;
samtools faidx fasta/all_species/cds_corrected_all_codons.fna \
-r id_lists/cds_palo_isoform_corrected/${orthogroup}.txt | \
seqkit replace -p '.*' -r '{nr}' > \
fasta/sco/cds_corrected/${orthogroup}.fna;
done

# for i in $(seq $(wc -l assembly_list.txt | cut -d ' ' -f 1)); do
# 	name=$(sort -k1 assembly_list.txt | sed -n ${i}p | cut -f 4)
# 	echo $name
# 	sed -i 's/>'${i}'$/>'${name}'/g' fasta/sco/cds_corrected/*.fna
# done

exit
