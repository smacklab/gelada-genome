#!/bin/bash

i=$1
j=$2

source scripts/_include_options.sh

int=$(printf %04d $(grep -n $j data/${genome}_regions.txt | cut -d ':' -f 1))
chr=$(printf %02d $(grep -n $(echo $j | cut -d ':' -f 1) ${genome_path}.fai | cut -d ':' -f 1))

mkdir -p bam

samtools view -b bam/${i}.${genome}.bam ${j} > bam/${i}.${genome}.region.${int}.chr${chr}.bam
samtools index -b bam/${i}.${genome}.region.${int}.chr${chr}.bam
