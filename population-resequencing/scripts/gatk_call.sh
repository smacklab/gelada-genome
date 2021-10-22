#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="gatk-hap"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/samples.txt)
s=$1

source scripts/_include_options.sh

if [ ! -z $2 ] && [ $2 -gt 0 ]; then
regions=$(echo $(tac data/${genome}_regions.txt))
else
regions=$(echo $(cat data/${genome}_regions.txt))
fi

module load perl/5.26.0
module load parallel/20180922
module load gatk/4.1.2.0
module load java/8u92
module load htslib/1.9.0

parallel -j $SLURM_CPUS_ON_NODE scripts/gatk_call_variants_single.sh {1} {2} {3} ::: $id ::: $regions ::: $s

exit

