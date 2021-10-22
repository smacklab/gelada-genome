#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="gatk-gvc"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

source scripts/_include_options.sh

int=$(printf %04d $SLURM_ARRAY_TASK_ID)
region=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/${genome}_regions.txt)
chr=$(printf %02d $(grep -n $(echo $region | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))
s=$1

module load perl/5.26.0
module load gatk/4.1.2.0
module load java/8u92
module load htslib/1.9.0

skip=0

if [ ! -f vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz ]; then
mkdir -p db
rm -rf db/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.step${s}
gatk --java-options "-Xmx100g" \
    GenomicsDBImport \
    --intervals ${region} \
    $(ls gvcf/*.${genome}.region.${int}.chr${chr}.raw.step${s}.g.vcf.gz | sed 's/^/--variant /g') \
    --batch-size 50 \
    --tmp-dir=/tmp \
    --reader-threads $SLURM_CPUS_ON_NODE \
    --genomicsdb-workspace-path db/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.step${s}

touch vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz

skip=1
fi

if [ $skip -lt 1 ]; then
if [ ! -f vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz.tbi ]; then
mkdir -p vcf

tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)

gatk --java-options "-Xmx100g" \
    GenotypeGVCFs \
    --intervals ${region} \
    --reference $genome_path \
    --variant gendb://db/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.step${s} \
    --output ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz

if [ $(bgzip -t ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz && echo 0 || echo 1) -eq 0 ]; then
mv ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz
mv ${tmp_dir}/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz.tbi vcf/${dataset}.${genome}.bootstrap.region.${int}.chr${chr}.raw.step${s}.vcf.gz.tbi
fi
fi
fi

exit
