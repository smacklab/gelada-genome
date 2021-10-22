#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="maf"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=1:00:00
#SBATCH --mem=2GB

cmd_prefix='wget http://ftp.ensembl.org/pub/release-101/maf/ensembl-compara/multiple_alignments/57_mammals.epo/'
maf=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/epo-maf.txt)

cmd=${cmd_prefix}${maf}

cd maf

$cmd

exit