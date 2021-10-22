#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=kchiou@uw.edu
#SBATCH --job-name="palo"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --account=csde
#SBATCH --partition=csde
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8GB

module load contrib/parallel/20171122

srun="srun --exclusive -N1 -n1"

mkdir -p logs

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog logs/parallel_${SLURM_JOB_ID}.log --resume"

$parallel "$srun scripts/run_palo.sh {1} > logs/parallel.${SLURM_JOB_ID}_{1}.log" ::: {4105..6807}

exit