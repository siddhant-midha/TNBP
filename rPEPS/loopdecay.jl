using JLD2
include("randompeps.jl")

# Parse command line arguments for SLURM job array
if length(ARGS) < 1
    println("Usage: julia loopdecay.jl <task_id>")
    exit(1)
end

task_id = parse(Int, ARGS[1])

# Define parameter combinations
N_values = [5, 10, 20]
eta_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

# Calculate total combinations
total_N = length(N_values)
total_eta = length(eta_values)

# Convert task_id to 0-based indexing
task_idx = task_id - 1

# Calculate which (N, eta) combination this task should run
N_idx = task_idx % total_N + 1  # Julia uses 1-based indexing
eta_idx = task_idx ÷ total_N + 1

# Get actual values
N = N_values[N_idx]
η = eta_values[eta_idx]

println("Task ID: $task_id -> N=$N, η=$η (N_idx=$N_idx, eta_idx=$eta_idx)")

# Configuration
T = N 
nsamples = 10000  # Reduced for reasonable runtime
wmax = 8  # Use highest weight that contains all loops
ti = true 
orthog = false 
annealing = 0.75 
maxiter = 1000 

# Load cluster data to get loop structures
cluster_data = load_latest_cluster_file(N, wmax; bc = "open", save_dir = "../saved_clusters")
loop_objects = cluster_data.all_loops  # This contains loops of all weights up to wmax
all_loops = [loop_object.edges for loop_object in loop_objects]

println("Loaded $(length(all_loops)) loops from cluster data")
# Storage for results across different loop weights
results_by_weight = Dict{Int, Vector{Float64}}()

println("Processing η = $η with $nsamples samples")

# Average over multiple samples
for sample in 1:nsamples
    if sample % 1000 == 0
        println("  Sample $sample/$nsamples")
    end
    
    # Get random PEPS tensor network for this η
    tensors, _ = peps_controllable(N, T; η=η, ti=ti, type="complex")
    adj_mat, edges, links = BP.get_adj_mat(tensors)
    messages = BP.get_messages(tensors, edges, links; random_part = 0.01) 
    messages = BP.message_passing(tensors, messages, edges, adj_mat; α=annealing, max_iters=maxiter, diagnose=false, normalise=true)
    Z_list = BP.get_fixed_point_list(tensors, messages, adj_mat)
    tensors = BP.normalize_tensors(tensors, Z_list)
    
    # Calculate loop contributions and organize by actual loop weight (length)
    for loop in all_loops
        if all([e in edges for e in loop])
            loop_weight = length(loop)  # Actual weight is the loop length
            loop_contri = scalar(BP.loop_contribution(loop, messages, tensors, edges, links, adj_mat))
            
            if !haskey(results_by_weight, loop_weight)
                results_by_weight[loop_weight] = Float64[]
            end
            
            push!(results_by_weight[loop_weight], abs(loop_contri))  # Take absolute value
        end
    end
end

# Calculate average contributions for each weight
avg_results = Dict()
for (w, contributions) in results_by_weight
    if !isempty(contributions)
        avg_results[w] = mean(contributions)
        println("Weight $w: $(length(contributions)) loops, avg contribution = $(avg_results[w])")
    end
end

# Save results

results = Dict(
    "N" => N,
    "eta" => η,
    "nsamples" => nsamples,
    "wmax" => wmax,
    "avg_contributions_by_weight" => avg_results,
    "raw_contributions_by_weight" => results_by_weight,
    "total_loops_processed" => sum(length(v) for v in values(results_by_weight))
)

# Create results directory if it doesn't exist
results_dir = "loop_decay_results_complex"
if !isdir(results_dir)
    mkdir(results_dir)
end

filename = joinpath(results_dir, "loop_decay_N$(N)_eta$(replace(string(η), "." => "p")).jld2")
save(filename, "results", results)
println("Saved results to: $filename")