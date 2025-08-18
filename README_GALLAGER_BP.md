# Gallager BP Loop Correction Analysis

This directory contains code for running Gallager LDPC code simulations with and without loop corrections using SLURM job arrays.

## Files

- `gallagerBP.jl` - Main simulation script that runs both loop and no-loop decoders
- `gallager.slurm` - SLURM job array script for parallel execution
- `collect_gallager_results.jl` - Script to collect and analyze results after jobs complete
- `collect_results.slurm` - SLURM script to automatically run collection after main jobs finish
- `ldpc_tanner_loops.jl` - Tanner graph loop finding functions
- `functions/BP.jl` - Belief propagation implementation

## Usage

### 1. Submit the main job array

```bash
sbatch gallager.slurm
```

This will submit 20 jobs (task IDs 1-20), each computing error rates for a different physical error rate `p` value from 0.1 to 0.2.

### 2. Monitor job progress

```bash
squeue -u $USER
```

### 3. Collect results after jobs complete

#### Option A: Automatic collection
Edit `collect_results.slurm` to replace `JOBID` with the actual job array ID from step 1, then:

```bash
sbatch collect_results.slurm
```

#### Option B: Manual collection
After all jobs complete:

```bash
julia collect_gallager_results.jl
```

## Output Structure

### Individual Job Results
Each job creates files in `gallagerBP0results/` with naming pattern:
```
results_task{task_id}_n{n}_p{p_value}_{timestamp}.jld2
```

Each file contains:
- Simulation parameters (n, p, d_v, d_c, etc.)
- Results with loop corrections: `mean_error_rate_loops`, `std_error_rate_loops`
- Results without loops: `mean_error_rate_no_loops`, `std_error_rate_no_loops`

### Collected Results
The collection script creates:
- `collected_results_{timestamp}.jld2` - Organized data for all n and p values
- `error_rates_vs_p_{timestamp}.png/pdf` - Error rate plots
- `loop_improvement_{timestamp}.png/pdf` - Loop correction improvement plots

## Parameters

Current configuration:
- **Code dimensions**: d_v=3, d_c=4 (regular LDPC)
- **Code sizes**: n ∈ [48, 96]
- **Physical error rates**: p ∈ [0.1, 0.2] (20 points)
- **Statistics**: 50 batches × 50 samples = 2500 samples per (n,p) point
- **Loop order**: Up to 6th order loops

## Modifying Parameters

To change simulation parameters, edit the parameter section in `gallagerBP.jl`:

```julia
d_v = 3 
d_c = 4
num_batches = 50
samples_per_batch = 50
max_loop_order = 6
n_ar = [48, 96]
ps = LinRange(0.1, 0.2, 20)
```

**Important**: If you change the number of p values, update the SLURM array range in `gallager.slurm`:
```bash
#SBATCH --array=1-20  # Change 20 to match length(ps)
```

## Troubleshooting

### Check job status
```bash
# View running/pending jobs
squeue -u $USER

# Check specific job output
cat logs/eta_JOBID_TASKID.out
cat logs/eta_JOBID_TASKID.err
```

### Check results
```bash
# List result files
ls gallagerBP0results/results_task*

# Check a specific result file
julia -e "using JLD2; data = load(\"gallagerBP0results/results_task1_n48_p0.1_*.jld2\"); println(keys(data))"
```

### Rerun failed jobs
If some jobs fail, you can rerun specific task IDs:
```bash
sbatch --array=5,7,12 gallager.slurm  # Rerun only tasks 5, 7, and 12
```

## Expected Runtime

- Each task processes 2 code sizes × 2500 samples ≈ 5000 simulations
- Estimated runtime: 2-4 hours per task depending on system load
- Total walltime for all 20 tasks: ~8 hours (with 4 CPU cores per task)

## Output Analysis

The collected results allow analysis of:
1. **Threshold behavior**: Error rate vs physical error rate curves
2. **Loop correction efficacy**: Comparison between decoders with/without loops
3. **Finite-size effects**: Performance differences between n=48 and n=96
4. **Statistical significance**: Error bars from batch-to-batch variation
