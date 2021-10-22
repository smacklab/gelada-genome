#!/bin/bash
module load samtools/1.9
module load bedtools2/2.24.0
module load bcftools/1.10.2

SMP=`sed -n ${SLURM_ARRAY_TASK_ID}p ../samples`

aln=/scratch/klchiou/theropithecus_wgs/bam/
genome=panubis1
pos_mask=/scratch/nsnyderm/thero_wgs/msmc2/pos_mask_beds

## calculate mean depth (just chromosome 20)
meandepth=`samtools depth -r NC_044995.1 ${aln}/${SMP}.${genome}.bam |  awk '{sum+=$3} END { printf "%.0f", sum/NR}'`

# 50% min
minDP=`echo $(( ${meandepth} / 2))`  ## 50% of mean depth
# 250% max
maxDP=`echo $(( ${meandepth} / 2 * 5))`  ## 300% of mean depth

echo $meandepth
echo $minDP
echo $maxDP

module load anaconda/py3
source activate mosdepth

export MOSDEPTH_Q0=NO_COVERAGE  
export MOSDEPTH_Q1=LOW_COVERAGE  
export MOSDEPTH_Q2=CALLABLE      
export MOSDEPTH_Q3=HIGH_COVERAGE 

mosdepth -n --fast-mode -t 20 --quantize 0:1:${minDP}:${maxDP}: mosdepth_out/${SMP} ${aln}/${SMP}.${genome}.bam

source deactivate 

zgrep CALLABLE mosdepth_out/${SMP}.quantized.bed.gz > ${pos_mask}/${SMP}.posmask.bed

cat panubis1_chroms | parallel --verbose -j 20 "grep {} ${pos_mask}/${SMP}.posmask.bed > ${pos_mask}/${SMP}.{}.posmask.bed"

aln=/scratch/klchiou/theropithecus_wgs/bam/
genome=panubis1
pos_mask=/scratch/nsnyderm/thero_wgs/msmc2/pos_mask_beds

BAM_DIR=/scratch/klchiou/theropithecus_wgs/bam/
VCF_DIR=/scratch/nsnyderm/thero_wgs/msmc2/vcfs # path to /bams dir 
OUT_DIR=/scratch/nsnyderm/thero_wgs/msmc2/msmc_prep

SMP=`sed -n ${SLURM_ARRAY_TASK_ID}p ../samples`

## these steps remove those sites that are missing or don't pass filter in a specific sample
## and remove them from the pos_mask file (so these sites weren't covered enough in this sample)

cat panubis1_chroms | parallel --verbose -j 20 "bcftools view -m2 -M2 -v snps -f .,PASS -s ${SMP} ${VCF_DIR}/gelada.{}.vcf.gz| \
   bcftools view -i 'F_MISSING<0.1' -Oz -o ${VCF_DIR}/${SMP}.{}.vcf.gz  -
subtractBed -a ${VCF_DIR}/gelada.{}.vcf.gz -b ${VCF_DIR}/${SMP}.{}.vcf.gz > nopass.${SMP}.{}.vcf
bcftools view -h ${VCF_DIR}/gelada.{}.vcf.gz > ${SMP}.{}.header
cat ${SMP}.{}.header nopass.${SMP}.{}.vcf > ${SMP}.{}.vcf
mv ${SMP}.{}.vcf nopass.${SMP}.{}.vcf 
subtractBed -a ${pos_mask}/${SMP}.{}.posmask.bed -b nopass.${SMP}.{}.vcf > ${SMP}.{}.posmask.bed
rm nopass.${SMP}.{}.vcf ${SMP}.{}.vcf ${SMP}.{}.header
mv ${SMP}.{}.posmask.bed ${pos_mask}/${SMP}.{}.posmask.bed"