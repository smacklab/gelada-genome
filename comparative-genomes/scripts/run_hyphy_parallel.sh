#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=chiou@asu.edu
#SBATCH --job-name="hyphy"
#SBATCH --output=out/slurm-%j.out
#SBATCH --error=out/slurm-%j.err
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB

# module load contrib/parallel/20171122

module load perl/5.26.0
module load parallel/20180922

srun="srun --exclusive -N1 -n1"

start=$(( ($SLURM_ARRAY_TASK_ID - 1) * 100 + 1 ))
end=$(( $start + 99 ))

if [ $end -gt $(( $(wc -l id_lists/orthogroups/sco_most.txt | cut -d ' ' -f 1) * 2 )) ]; then
	end=$(( $(wc -l id_lists/orthogroups/sco_most.txt | cut -d ' ' -f 1) * 2 ))
fi

mkdir -p logs

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog logs/parallel_${SLURM_JOB_ID}.log --resume"

$parallel "$srun scripts/run_hyphy_single.sh {1} > logs/parallel.${SLURM_JOB_ID}_{1}.log" ::: $(eval echo -e {${start}..${end}})

exit