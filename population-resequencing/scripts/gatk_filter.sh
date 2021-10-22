#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="gatk-vfl"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

s=$1

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

source scripts/_include_options.sh

module load perl/5.26.0
module load parallel/20180922
module load gatk/4.1.2.0
module load java/8u92
module load bcftools/1.9
module load tabix/0.2.6
module load htslib/1.9.0

parallel -j $slots scripts/gatk_filter_variants_single.sh {1} {2} ::: $(eval echo -e {1..$(wc -l data/${genome}_regions.txt | cut -d ' ' -f 1)}) ::: $s

exit
