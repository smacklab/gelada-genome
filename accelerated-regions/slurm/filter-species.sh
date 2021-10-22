#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="spefilter"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=4:00:00
#SBATCH --exclusive

module load perl/5.26.0
module load parallel/20180922

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

rm -rf data/maf-pass.txt
rm -rf data/maf-fail.txt

parallel -j $slots scripts/filter-species.sh {1} ::: $(eval echo -e {1..$(wc -l data/maf-files.txt | xargs | cut -d ' ' -f 1)})

exit