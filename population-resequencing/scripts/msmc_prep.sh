#!/bin/bash
module load msmc2/2019.2 
module load python/3.7.1

pos_mask=/scratch/nsnyderm/thero_wgs/msmc2/pos_mask_beds
VCF_DIR=/scratch/nsnyderm/thero_wgs/msmc2/vcfs
OUT_DIR=/scratch/nsnyderm/thero_wgs/msmc2/msmc_prep
MSMC_TOOLS=/scratch/nsnyderm/thero_wgs/msmc2/msmc-tools                         

SMP=`sed -n ${SLURM_ARRAY_TASK_ID}p ../samples`

cat panubis1_chroms | parallel --verbose -j 10 "${MSMC_TOOLS}/generate_multihetsep.py ${VCF_DIR}/${SMP}.{}.vcf.gz \
--mask ${pos_mask}/${SMP}.{}.posmask.bed > ${OUT_DIR}/${SMP}.{}.postMultiHetSep.txt"