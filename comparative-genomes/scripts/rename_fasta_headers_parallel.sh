#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="rehead"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB

# module load contrib/parallel/20171122
module load parallel/20180922

srun="srun --exclusive -N1 -n1"

mkdir -p logs

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog logs/parallel_${SLURM_JOB_ID}.log --resume"

$parallel "$srun scripts/rename_fasta_headers.sh {1} > logs/parallel.${SLURM_JOB_ID}_{1}.log" ::: {1..7330}

exit