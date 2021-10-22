#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="vcfconcat"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

s=$1

source scripts/_include_options.sh

module load perl/5.26.0
module load parallel/20180922
module load bcftools/1.9

if [ $(cut -d ':' -f 1 data/${genome}_regions.txt | sort -u | wc -l | xargs) -gt $SLURM_CPUS_ON_NODE ]; then
slots=$SLURM_CPUS_ON_NODE
else
slots=$(cut -d ':' -f 1 data/${genome}_regions.txt | sort -u | wc -l | xargs)
fi

parallel -j $slots scripts/bcftools_concat_variants_single.sh {1} {2} ::: $(eval echo -e {01..$(cut -d ':' -f 1 data/${genome}_regions.txt | sort -u | wc -l | xargs)}) ::: $s

exit
