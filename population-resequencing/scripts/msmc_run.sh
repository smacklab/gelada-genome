#!/bin/bash

module load msmc2/2019.2 

INPUT_DIR=/scratch/nsnyderm/thero_wgs/msmc2/msmc_prep
OUT_DIR=/scratch/nsnyderm/thero_wgs/msmc2/msmc_final
MSMC_TOOLS=/scratch/nsnyderm/thero_wgs/msmc2/msmc-tools                         

SMP=`sed -n ${SLURM_ARRAY_TASK_ID}p ../samples`

msmc2 -t 8 -o ${OUT_DIR}/${SMP}.msmc ${INPUT_DIR}/${SMP}*postMultiHetSep.txt