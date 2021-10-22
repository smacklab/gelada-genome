#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="treefam"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=2-00:00:00
#SBATCH --mem=48GB

module load perl/5.26.0
# Do not use hmmer version > 3.0
export PATH=/data/nsnyderm/software/hmmer/3.0/bin:$PATH

export PERL5LIB=$(pwd)/treefam_tools/treefam_scan:/home/klchiou/perl5/lib/perl5:$PERL5LIB

accession=$(sed -n ${SLURM_ARRAY_TASK_ID}p assembly_list.txt | cut -f 1 | sed 's/ /_/g')
genome=$(sed -n ${SLURM_ARRAY_TASK_ID}p assembly_list.txt | cut -f 2 | sed 's/ /_/g')

cd treefam_tools/treefam_scan

fasta=fasta/${accession}_${genome}_protein.faa

mkdir -p results

perl treefam_scan.pl -fasta $fasta -dir hmm_lib/ -hmm_file TreeFam9 -cpu $SLURM_CPUS_ON_NODE -outfile results/treefam_hmm_${genome}_results.txt

exit;`