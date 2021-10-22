#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="hierfstats"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --partition=mrline,mrline-serial,gmascaro,cidsegpu3,epyc1,spressecpu1,lsankargpu1,lsankargpu2,gdcsgpu1,rcgpu3
#SBATCH --qos=wildfire

p=$1 # prefix
s=$2 # step

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

source scripts/_include_options.sh

module load perl/5.26.0
module load parallel/20180922
module load r/4.0.2

# vcf/gelada.tgel1.bootstrap.region.0001.chr01.pas.step2.vcf.gz

parallel -j $slots scripts/hierfstat_stats.R {1} {2} {3} ::: $p ::: $(eval echo -e {1..$(wc -l data/${genome}_regions.txt | cut -d ' ' -f 1)}) ::: $s

# vcf/gelada.tgel1.bootstrap.region.0001.chr01.pas.step2.fst.rds
# scripts/hierfstat_stats.sh gelada.tgel1 2

# sbatch --nodes=2 scripts/hierfstat_stats.sh gelada.tgel1 2