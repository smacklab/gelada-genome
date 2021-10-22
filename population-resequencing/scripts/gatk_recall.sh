#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="gatk-fxh"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

s=$1

ids=$(echo $(tac data/samples.txt))
r=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/regions.txt)

int=$(printf %04d $SLURM_ARRAY_TASK_ID)

rm -rf gvcf/*.panubis1.region.${int}.chr*.raw.step${s}.g.vcf.gz
rm -rf gvcf/*.panubis1.region.${int}.chr*.raw.step${s}.g.vcf.gz.tbi

module load perl/5.26.0
module load parallel/20180922
module load gatk/4.1.2.0
module load java/8u92
module load htslib/1.9.0

parallel -j $SLURM_CPUS_ON_NODE scripts/gatk_call_variants_single.sh {1} {2} {3} ::: $ids ::: $r ::: $s

exit
