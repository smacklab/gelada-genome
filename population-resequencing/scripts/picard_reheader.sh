#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="gatk"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/samples.txt)

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

source scripts/_include_options.sh

module load picard/2.9.2
module load perl/5.26.0
module load parallel/20180922
module load gatk/4.1.2.0
# module load java/8u92
module load samtools/1.9

parallel -j $slots scripts/picard_reheader_single.sh {1} {2} ::: $id ::: $(echo $(cat data/${genome}_regions.txt))

exit

