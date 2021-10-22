#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="mitocaller"
#SBATCH --output=out/slurm_%A-%a.out
#SBATCH --error=out/slurm_%A-%a.err
#SBATCH --time=7-00:00:00
#SBATCH --partition=mrline,mrline-serial,gmascaro,cidsegpu3,epyc1,spressecpu1,lsankargpu1,lsankargpu2,gdcsgpu1,rcgpu3
#SBATCH --qos=wildfire

export PATH=$PATH:/scratch/klchiou/sw/mitocaller/bin

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/samples.txt)

module load samtools/1.9

source scripts/_include_options.sh

chrmt=$(tail -n1 $genome_path.fai | cut -f 1)

mkdir -p mt-bam

samtools view -bh bam/${id}.${genome}.bam $chrmt > mt-bam/${id}.mt.bam

mkdir -p mt-calls

module purge
module load gcc/6.3.0
mitoCaller -m -b mt-bam/${id}.mt.bam -r $genome_path > mt-calls/${id}.summary.txt

module purge

module load r/4.0.2

mkdir -p cytb

scripts/call_cytb.R $id > cytb/${id}.cytb.fa

# cytb = subset(foo,Pos >=14177 & Pos < 15910)
# identical(cytb$Pos,14177:15909)
# paste(apply(matrix(as.integer(gsub('[ACGT]:','',as.matrix(cytb[,c('FilteredNumA','FilteredNumC','FilteredNumG','FilteredNumT')]))),ncol=4),1,function(x) if ((x[which.max(x)] - max(x[setdiff(1:4,which.max(x))])) <= 0) 'N' else c('A','C','G','T')[which.max(x)]),collapse='')