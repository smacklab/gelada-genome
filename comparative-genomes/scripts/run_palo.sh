#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=kchiou@uw.edu
#SBATCH --job-name="palo"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --account=csde-ckpt
#SBATCH --partition=ckpt
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB

wd=$(pwd)

orthogroup=$(sed -n ${SLURM_ARRAY_TASK_ID}p id_lists/orthogroups/sco_all.txt)

module load contrib/palo/20130118
module load contrib/r/3.6.1

while [ $(wc -w ${wd}/palo-tmp/${orthogroup}/to_do.txt | cut -d ' ' -f 1) -gt 0 ]; do
# if [ $(wc -w ${wd}/palo-tmp/${orthogroup}/to_do.txt | cut -d ' ' -f 1) -gt 0 ]; then echo 'Go!'; fi
	cd ${wd}
	scripts/prep_palo_input_files.R ${orthogroup}
	cd palo-tmp/${orthogroup}
	palo2.py
	if [ ! -f files/palo_combs.txt ] || [ $(wc -l files/palo_combs.txt | cut -d ' ' -f 1) -lt 1 ]; then exit; fi
	mv to_do_next.txt to_do.txt
	if [ $(wc -l ${wd}/palo-tmp/${orthogroup}/files/unprocessed_genes.txt | cut -d ' ' -f 1) -gt 0 ]; then break; fi
done

cd $wd

scripts/summarize_palo.R ${orthogroup}

module load contrib/samtools/1.9

mkdir -p fasta/sco
mkdir -p fasta/sco/aa
mkdir -p fasta/sco/cds

module load contrib/seqkit/0.10.0

samtools faidx fasta/all_species/proteins.faa \
	-r id_lists/aa_palo_isoform/${orthogroup}.txt | \
	seqkit replace -p '.*' -r '{nr}' > \
	fasta/sco/aa/${orthogroup}.faa

samtools faidx fasta/all_species/cds.fna \
	-r id_lists/cds_palo_isoform/${orthogroup}.txt | \
	seqkit replace -p '.*' -r '{nr}' > \
	fasta/sco/cds/${orthogroup}.fna

exit
