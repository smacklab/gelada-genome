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


module load samtools/1.9

source scripts/_include_options.sh

chrmt=$(tail -n1 $genome_path.fai | cut -f 1)

mkdir -p mt-sam
mkdir -p mt-fq
mkdir -p mt-fq-singletons

samtools view -h bam/${id}.${genome}.bam $chrmt > mt-sam/${id}.mt.sam

samtools fastq mt-sam/${id}.mt.sam -s mt-fq-singletons/{id}.fastq -1 mt-fq/${id}.R1.fastq -2 mt-fq/${id}.R2.fastq

module load python/3.7.1
module load spades/3.15.2
module load bowtie2/2.4.1
module load ncbi-blast/2.6.0

get_organelle_from_reads.py -1 fastq/${id}.R1.fastq.gz -2 mt-fq/${id}.R2.fastq.gz -s genomes/FJ785426.1.fa -F animal_mt -o mt/${id} -R 10 --reduce-reads-for-coverage Inf # -t $SLURM_CPUS_ON_NODE