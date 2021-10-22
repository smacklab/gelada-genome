#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="gatk-bsr"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

# --ntasks=1
# --cpus-per-task=40
# --mem-per-cpu=4GB

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/samples.txt)
s=$1

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

source scripts/_include_options.sh

module load perl/5.26.0
module load parallel/20180922
module load gatk/4.1.2.0
module load java/8u92
module load samtools/1.9

parallel -j $slots scripts/gatk_recal_single.sh {1} {2} {3} ::: $id ::: $(echo $(cat data/${genome}_regions.txt)) ::: $s

exit

