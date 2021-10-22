#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="extract"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=2GB

module load perl/5.26.0
module load parallel/20180922

scripts/extract-blocks.sh $SLURM_ARRAY_TASK_ID $1 $2 $3

exit
