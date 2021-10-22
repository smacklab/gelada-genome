#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="gatk-gvc"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --partition=gpu
#SBATCH --qos=wildfire
#SBATCH --gres=gpu:4
#SBATCH --constraint=V100
#SBATCH --ntasks=32

source scripts/_include_options.sh

int=$(printf %04d $SLURM_ARRAY_TASK_ID)
region=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/${genome}_regions.txt)
chr=$(printf %02d $(grep -n $(echo $region | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))
s=$1

module load parabricks/3.5.0
module load htslib/1.9.0

mkdir -p gdb

pbrun importgvcftodb \
	--db-dir gdb/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.step${s} \
	$(/usr/bin/ls gvcf/*.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz | sed 's/^/--in-gvcf /g')

mkdir -p var/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.step${s}

pbrun selectvariants \
	--ref $genome_path \
	--db-dir gdb/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.step${s} \
	--out-gvcf-dir var/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.step${s}

tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)

pbrun genotypegvcf \
	--ref $genome_path \
	--in-selectvariants-dir var/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.step${s} \
	--out-vcf ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz

mkdir -p parabricks/vcf
if [ $(bgzip -t ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz && echo 0 || echo 1) -eq 0 ]; then
mv ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz parabricks/vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz
mv ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz.tbi parabricks/vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz.tbi
fi

exit
