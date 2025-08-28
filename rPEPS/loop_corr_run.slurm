#!/bin/bash
#SBATCH --job-name=LoopCorrRun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=logs/loop_corr_%A_%a.out
#SBATCH --error=logs/loop_corr_%A_%a.err
#SBATCH --array=1-120
#SBATCH --time=10:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=sm7456@princeton.edu

module purge
module load julia/1.11.1
export JULIA_PROJECT="/home/sm7456/.julia/itensorenv"
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Starting loop correction job array task $SLURM_ARRAY_TASK_ID"
echo "Job started at: $(date)"

julia loop_corr_run.jl $SLURM_ARRAY_TASK_ID

echo "Job completed at: $(date)"
