#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="download"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=2:00:00
#SBATCH --exclusive

module load perl/5.26.0
module load parallel/20180922

parallel -j $SLURM_CPUS_ON_NODE scripts/download-maf.sh {1} ::: $(eval echo -e {1..$(wc -l data/maf-files.txt | xargs | cut -d ' ' -f 1)})

exit