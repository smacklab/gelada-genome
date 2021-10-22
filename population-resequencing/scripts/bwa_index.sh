#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="bwa"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

module load bwa/0.7.17

source scripts/_include_options.sh

if [ -f ${genome_path}.gz ]; then gunzip ${genome_path}.gz; fi

bwa index $genome_path

module load samtools/1.9

samtools faidx $genome_path

module load gatk/4.1.2.0
module load java/8u92

gatk --java-options '-Xmx4g' CreateSequenceDictionary -R $genome_path