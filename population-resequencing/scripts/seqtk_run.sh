#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="seqtk"
#SBATCH --output=out/slurm_%A-%a.out
#SBATCH --error=out/slurm_%A-%a.err
#SBATCH --time=7-00:00:00
#SBATCH --partition=mrline,mrline-serial,gmascaro,cidsegpu3,epyc1,spressecpu1,lsankargpu1,lsankargpu2,gdcsgpu1,rcgpu3
#SBATCH --qos=wildfire

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/samples.txt)

mkdir -p fq

module load seqtk/1.3

seqtk sample -s100 -2 fastq/${id}.R${1}.fastq.gz 10000000 | gzip -c > fq/${id}.R${1}.fastq.gz
