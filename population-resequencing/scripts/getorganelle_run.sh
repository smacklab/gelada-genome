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
#SBATCH --exclusive

module load python/3.7.1
module load spades/3.15.2
module load bowtie2/2.4.1
module load ncbi-blast/2.6.0

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p data/samples.txt)

# mkdir -p fq

# module load seqtk/1.3

# seqtk sample -s100 -2 fastq/${id}.R${1}.fastq.gz 2000000 | gzip -c > fq/${id}.R${1}.fastq.gz

# The below line needs to be run first
# get_organelle_config.py -a animal_mt

mkdir -p mt

# get_organelle_from_reads.py -1 fq/${id}.R1.fastq.gz -2 fq/${id}.R2.fastq.gz -F animal_mt -s genomes/FJ785426.1.fa -o mt/${id} -R 10 -t $SLURM_CPUS_ON_NODE # --reduce-reads-for-coverage Inf -t $SLURM_CPUS_ON_NODE
get_organelle_from_reads.py -1 fq/${id}.R1.fastq.gz -2 fq/${id}.R2.fastq.gz -F animal_mt -s genomes/FJ785426.1.fa -o mt/${id} -R 20 -t $SLURM_CPUS_ON_NODE -w 30 --max-reads inf --reduce-reads-for-coverage inf # --reduce-reads-for-coverage Inf -t $SLURM_CPUS_ON_NODE

# module load emboss/6.6.0

# h01.fa

# mkdir -p mt-aln

# if [ $(grep -c '>' mt-fasta/${id}.mt.fa) -gt 1 ]; then
# for i in $(seq $(grep -c '>' mt-fasta/${id}.mt.fa)); do
# sed -n $(grep -n '>' mt-fasta/${id}.mt.fa | sed -n ${i}p | cut -d ':' -f 1),$(($(grep -n '>' mt-fasta/${id}.mt.fa | sed -n ${i}p | cut -d ':' -f 1)+1))p mt-fasta/${id}.mt.fa > mt-fasta/${id}.mt.scaffold${i}.fa
# water -asequence mt-fasta/${id}.mt.scaffold${i}.fa -bsequence h01.fa -snucleotide1 -snucleotide2 -scircular1 -aformat3 fasta -gapopen 10 -gapextend 0.5 -outfile mt-aln/${id}.${i}.aln.fa
# done
# else
# water -asequence mt-fasta/${id}.mt.fa -bsequence h01.fa -snucleotide1 -snucleotide2 -scircular1 -aformat3 fasta -gapopen 10 -gapextend 0.5 -outfile mt-aln/${id}.aln.fa
# fi

# Load fasta tools
# https://github.com/b-brankovics/fasta_tools

# export PATH=$PATH:/scratch/klchiou/downloads/fasta_tools/bin

# genomes/FJ785426.1.fa

