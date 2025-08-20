#!/bin/bash
#SBATCH --job-name=GallagerBPLoops # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1              # total number of tasks across all nodes
#SBATCH --cpus-per-task=4       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                # total memory (across all tasks)
#SBATCH --output=logs/eta_%A_%a.out
#SBATCH --error=logs/eta_%A_%a.err
#SBATCH --array=1-20
#SBATCH --time=8:00:00           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=sm7456@princeton.edu

module purge
module load julia/1.11.1
export JULIA_PROJECT="/home/sm7456/.julia/itensorenv"
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia gallagerBP.jl $SLURM_ARRAY_TASK_ID
