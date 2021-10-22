#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=kchiou@uw.edu
#SBATCH --job-name="ortho"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --account=csde
#SBATCH --partition=csde
#SBATCH --time=336:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --mem=248GB

module load contrib/orthofinder/2.3.3

orthofinder.py -t 28 -o $(pwd)/orthogroups -f fasta/aa_longest_isoform

exit
