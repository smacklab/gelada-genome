#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="preprocess"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=4:00:00
#SBATCH --exclusive

module load perl/5.26.0
module load parallel/20180922

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

parallel -j $slots scripts/preprocess-maf.sh {1} ::: $(eval echo -e {1..$(wc -l data/maf-pass.txt | xargs | cut -d ' ' -f 1)})

exit