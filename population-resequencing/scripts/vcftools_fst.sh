#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="vcf-fst"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

s=$1

# i=$1
# int=$(printf %04d $i)
# region=$(sed -n ${i}p data/${genome}_regions.txt)
# chr=$(printf %02d $(grep -n $(echo $region | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))
# s=$2

slots=$(echo $SLURM_JOB_CPUS_PER_NODE | sed 's/[()]//g' | sed 's/x/*/g' | sed 's/,/+/g' | bc)

source scripts/_include_options.sh

if [ ! -f data/pop1.txt ]; then
cat data/samples.txt | grep -v '^FIL' > data/pop1.txt
fi
if [ ! -f data/pop2.txt ]; then
cat data/samples.txt | grep '^FIL' > data/pop2.txt
fi

# VCFtools should already be loaded

mkdir -p results
mkdir -p results/fst

window_size=1000
step_size=500

# Windowed
vcftools \
	--gzvcf vcf_final/${dataset}.${genome}.bootstrap.all.step${s}.vcf.gz \
	--weir-fst-pop data/pop1.txt --weir-fst-pop data/pop2.txt \
	--fst-window-size $window_size --fst-window-step $step_size \
	--out results/fst/${dataset}.${genome}.bootstrap.all.step${s}

# Per site
vcftools \
	--gzvcf vcf_final/${dataset}.${genome}.bootstrap.all.step${s}.vcf.gz \
	--weir-fst-pop data/pop1.txt --weir-fst-pop data/pop2.txt \
	--out results/fst/${dataset}.${genome}.bootstrap.all.step${s}_persite
#	--fst-window-size $window_size --fst-window-step $step_size \

# 
# vcftools \
# 	--gzvcf vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.vcf.gz \
# 	--weir-fst-pop data/pop1.txt --weir-fst-pop data/pop2.txt \
# 	--fst-window-size $window_size --fst-window-step $step_size \
# 	--out results/fst/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.pas.step${s}.fst
# 