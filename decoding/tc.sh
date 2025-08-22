#!/bin/zsh
#SBATCH --job-name=tc
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=23:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL          # send email when job ends
#SBATCH --array=0-399
#SBATCH --mail-user=yz4281@princeton.edu

let "id=1000+$SLURM_ARRAY_TASK_ID+1"
let date=250809

echo "id: $id"
echo "date: $date"

# Create date directory if it doesn't exist
mkdir -p data/$date

# Run Julia with proper filename
julia tc.jl "$date/$id"