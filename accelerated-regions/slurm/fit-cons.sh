#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="cons"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=4:00:00
#SBATCH --exclusive

module load perl/5.26.0
module load parallel/20180922

source config/settings.cfg

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

parallel -j $slots scripts/fit-cons.sh {1} {2} {3} {4} ::: $(eval echo -e {1..$(wc -l data/${test_species}-blocks.txt | xargs | cut -d ' ' -f 1)}) ::: $1 ::: $2 ::: $3

exit