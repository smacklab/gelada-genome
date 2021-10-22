#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="preprocess"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=16GB

module load perl/5.26.0
module load parallel/20180922

scripts/preprocess-maf.sh $SLURM_ARRAY_TASK_ID

exit