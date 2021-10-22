#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="speedseq"
#SBATCH --output=out/slurm-%j_%a.out
#SBATCH --error=out/slurm-%j_%a.err
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

module load speedseq/0.1.2

source scripts/_include_options.sh

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/samples.txt)

speedseq align -v -t $SLURM_CPUS_ON_NODE -M 150 -R "@RG\tID:$id\tSM:$id\tLB:$id" \
	-o bam/${id}.${genome} $genome_path \
	fastq/${id}.R1.fastq.gz fastq/${id}.R2.fastq.gz
