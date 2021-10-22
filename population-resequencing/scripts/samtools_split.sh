#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="sam-view"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/samples.txt)

source scripts/_include_options.sh

module load perl/5.26.0
module load parallel/20180922
module load samtools/1.9

parallel -j $SLURM_CPUS_ON_NODE scripts/samtools_extract_region.sh {1} {2} ::: $id ::: $(echo $(cat data/${genome}_regions.txt))

exit

