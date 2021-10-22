# Mapping and genotyping pipeline

Scripts in this folder cover analyses performed on population resequencing data.

The steps below outline the code run to perform mapping and genotyping from population resequencing data. It can be run and repeated with multiple configurations by editing the `scripts/_include_options.sh` file.

```
sbatch --array=1-90 scripts/speedseq-align.sh
sbatch --array=1-90 scripts/samtools_split.sh
sbatch --array=1-90 scripts/picard_reheader.sh

sbatch --array=1-90 scripts/gatk_call.sh 0
sbatch --array=1-$(wc -l data/${genome}_regions.txt | cut -d ' ' -f 1) scripts/gatk_genotype.sh 0
sbatch scripts/gatk_filter.sh 0

sbatch --array=1-90 scripts/gatk_recal.sh 1
sbatch --array=1-90 scripts/gatk_call.sh 1
sbatch --array=1-$(wc -l data/${genome}_regions.txt | cut -d ' ' -f 1) scripts/gatk_genotype.sh 1
sbatch scripts/gatk_filter.sh 1
                             
sbatch --array=1-90 scripts/gatk_recal.sh 2
sbatch --array=1-90 scripts/gatk_call.sh 2 
sbatch --array=1-$(wc -l data/${genome}_regions.txt | cut -d ' ' -f 1) scripts/gatk_genotype.sh 2     
sbatch scripts/gatk_filter.sh 2

# Concatenate VCFs and assemble all the data
sbatch scripts/bcftools_concat.sh 2

sbatch scripts/vcftools_filter.sh 2
```