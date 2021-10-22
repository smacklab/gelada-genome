#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=klchiou@asu.edu
#SBATCH --job-name="clustalo"
#SBATCH --output=out/slurm_%A-%a.out
#SBATCH --error=out/slurm_%A-%a.err
#SBATCH --time=7-00:00:00
#SBATCH --partition=mrline,mrline-serial,gmascaro,cidsegpu3,epyc1,spressecpu1,lsankargpu1,lsankargpu2,gdcsgpu1,rcgpu3
#SBATCH --qos=wildfire
#SBATCH --exclusive

export PATH=$PATH:/scratch/klchiou/sw/clustalomega/1.2.4/bin/

cat cytb/*.cytb.fa zinner.fa > cytb_all.fa

clustalo --threads=$SLURM_CPUS_ON_NODE --force -i cytb_all.fa -o cytb_all.aln.fa --guidetree-out cytb_all.nwk
clustalo --threads=$SLURM_CPUS_ON_NODE --force -i cytb_haplotypes.fa -o cytb_haplotypes.aln.fa --guidetree-out cytb_haplotypes.nwk

clustalo --threads=$SLURM_CPUS_ON_NODE --force -i cytb_all.fa --outfmt=phy -o cytb_all.aln.phy --guidetree-out cytb_all.nwk
clustalo --threads=$SLURM_CPUS_ON_NODE --force -i cytb_haplotypes.fa --outfmt=phy -o cytb_haplotypes.aln.phy --guidetree-out cytb_haplotypes.nwk