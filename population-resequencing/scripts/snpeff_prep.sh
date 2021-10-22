#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="snpeff-prep"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

source scripts/_include_options.sh

mkdir -p snpEff
mkdir -p snpEff/data
mkdir -p snpEff/data/genomes
mkdir -p snpEff/data/${genome}

if [ ! -f snpEff/data/genomes/${genome}.fa ]; then
echo 'Copying genome'
cp $genome_path snpEff/data/genomes/${genome}.fa
fi

if [ ! -f snpEff/gelada.bootstrap.all.step2.vcf.gz ]; then
echo 'Copying VCF file'
cp vcf_all/gelada.bootstrap.all.step2.vcf.gz snpEff/gelada.bootstrap.all.step2.vcf.gz
fi

if [ ! -f snpEff/data/panubis1/genes.gtf ]; then
echo 'Downloading and unpacking gtf'
wget -O snpEff/data/panubis1/genes.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/728/515/GCF_008728515.1_Panubis1.0/GCF_008728515.1_Panubis1.0_genomic.gtf.gz
gunzip snpEff/data/panubis1/genes.gtf.gz
fi

cd snpEff
echo "panubis1.genome : panubis1" > snpEff.config

module load snpEff/4.3

echo 'Building'
snpEff build -c ./snpEff.config -dataDir ./data -gtf22 -v panubis1

module load htslib/1.9.0

echo 'Annotating'
snpEff ann -c snpEff.config -dataDir data panubis1 gelada.bootstrap.all.step2.vcf.gz | bgzip -@ $SLURM_CPUS_ON_NODE > gelada.bootstrap.all.step2.ann.vcf.gz

echo 'Indexing'
tabix gelada.bootstrap.all.step2.ann.vcf.gz