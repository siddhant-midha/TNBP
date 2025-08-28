include("randompeps.jl")
include("cluster_exp.jl")
using JLD2, FileIO, Dates, Statistics


function run_cluster_correction_analysis(N::Int, w::Int, η::Float64, nsamples::Int)
    """
    Run cluster correction analysis for given parameters.
    Uses complex arithmetic in log calculations to avoid domain errors.
    Computes free energy per site errors: log(Z) / (N*T).
    
    Now loads unique_global_clusters directly from saved files when available,
    avoiding the need for manual cluster deduplication during analysis.
    Falls back to clusters_by_site with manual deduplication for older saved files.
    
    Returns Dict with results including average free energy per site errors.
    """
    
    println("="^80)
    println("🔥 Running Cluster Correction Analysis")
    println("Parameters: N=$N, w=$w, η=$η, nsamples=$nsamples")
    println("="^80)
    
    # Load cluster data once
    cluster_data = nothing
    unique_global_clusters = nothing
    clusters_by_site = nothing  # Fallback for old format
    all_loops = nothing
    
    try
        cluster_data = load_latest_cluster_file(N, w; bc = "open", save_dir = "../saved_clusters")
        # Try to load unique global clusters from the new format        
        println("✅ Loaded cluster data successfully")
    catch e
        println("❌ Failed to load cluster data for N=$N, w=$w: $e")
        return nothing
    end
    
    # Storage for results
    bp_errors = Float64[]           # Free energy per site errors
    cluster_errors = Float64[]      # Free energy per site errors with cluster correction
    bp_fe_densities = Float64[]     # BP free energy per site values
    exact_fe_densities = Float64[]  # Exact free energy per site values
    cluster_corrections = Float64[] # FE density corrections (real parts)
    
    # Fixed parameters
    ti = true   
    orthog = false   
    normalise = true 
    annealing = 0.75 
    maxiter = 1000
    T = N
    
    println("🔄 Running $nsamples samples...")
    
    for sample in 1:nsamples
        if sample % 50 == 0 || sample <= 10
            println("  Sample $sample/$nsamples")
        end
        
        try
            # Generate random PEPS
            tensors, peps = peps_controllable(N, T; η=η, ti=ti, type="complex")
            exact_PF = contract_peps_no_phys(peps; cutoff=1E-12, maxdim=2^N)
            exact_FE_per_site = log(Complex(exact_PF)) / (N*T)  # Complex log to handle any edge cases
            
            # BP calculation - ensure complex arithmetic throughout
            adj_mat, edges, links = BP.get_adj_mat(tensors)
            messages = BP.get_messages(tensors, edges, links; random_part=0.01)
            messages = BP.message_passing(tensors, messages, edges, adj_mat; α=annealing, max_iters=maxiter)
            
            Z_l = BP.get_fixed_point_list(tensors, messages, adj_mat)
            # Ensure Z_l entries are complex before taking log
            Z_l_complex = Complex.(Z_l)
            bp_FE_per_site = sum(log.(Z_l_complex)) / (N*T)
            
            T_normalized = BP.normalize_tensors(tensors, Z_l)
            
            # Cluster correction using unique global clusters
            clustercorrx = cluster_contr_by_site(cluster_data, T_normalized, messages, edges, links, adj_mat)
            
            cluster_FE_per_site_correction = clustercorrx / (N*T)
            corrected_FE_per_site = bp_FE_per_site + cluster_FE_per_site_correction
            
            # Calculate free energy per site errors
            bp_error = abs(bp_FE_per_site - exact_FE_per_site)
            cluster_error = abs(corrected_FE_per_site - exact_FE_per_site)
            
            # Store results
            push!(bp_errors, bp_error)
            push!(cluster_errors, cluster_error)
            push!(bp_fe_densities, real(bp_FE_per_site))
            push!(exact_fe_densities, real(exact_FE_per_site))
            push!(cluster_corrections, real(cluster_FE_per_site_correction))
            
        catch e
            println("⚠️  Error in sample $sample: $e")
            continue
        end
    end
    
    if isempty(bp_errors)
        println("❌ No successful samples!")
        return nothing
    end
    
    # Calculate statistics
    results = Dict(
        "N" => N,
        "w" => w,
        "eta" => η,
        "nsamples_requested" => nsamples,
        "nsamples_successful" => length(bp_errors),
        "bp_fe_error_mean" => mean(bp_errors),
        "bp_fe_error_std" => std(bp_errors),
        "cluster_fe_error_mean" => mean(cluster_errors),
        "cluster_fe_error_std" => std(cluster_errors),
        "fe_improvement_mean" => mean(bp_errors) - mean(cluster_errors),
        "fe_improvement_std" => std(bp_errors .- cluster_errors),
        "cluster_correction_mean" => mean(cluster_corrections),
        "cluster_correction_std" => std(cluster_corrections),
        "timestamp" => string(now()),
        # Store all raw data for detailed analysis
        "bp_fe_errors" => bp_errors,
        "cluster_fe_errors" => cluster_errors,
        # Keep old names for backward compatibility
        "bp_errors" => bp_errors,
        "cluster_errors" => cluster_errors,
        "bp_fe_densities" => bp_fe_densities,
        "exact_fe_densities" => exact_fe_densities,
        "cluster_corrections" => cluster_corrections
    )
    
    println("✅ Analysis complete!")
    println("  Successful samples: $(length(bp_errors))/$nsamples")
    println("  BP free energy per site error (mean ± std): $(mean(bp_errors)) ± $(std(bp_errors))")
    println("  Cluster free energy per site error (mean ± std): $(mean(cluster_errors)) ± $(std(cluster_errors))")
    println("  Free energy per site improvement: $(mean(bp_errors) - mean(cluster_errors))")
    
    return results
end

function get_parameter_combination(task_id::Int, N_list::Vector{Int}, w_list::Vector{Int}, η_list::Vector{Float64})
    """
    Convert task_id to (N, w, η) combination.
    Task IDs start from 1.
    """
    total_combinations = length(N_list) * length(w_list) * length(η_list)
    
    if task_id < 1 || task_id > total_combinations
        error("Task ID $task_id is out of range [1, $total_combinations]")
    end
    
    # Convert to 0-based indexing for easier calculation
    idx = task_id - 1
    
    # Unravel the index
    n_η = length(η_list)
    n_w = length(w_list)
    
    η_idx = idx % n_η + 1
    w_idx = (idx ÷ n_η) % n_w + 1  
    N_idx = idx ÷ (n_η * n_w) + 1
    
    return N_list[N_idx], w_list[w_idx], η_list[η_idx]
end

function save_results(results::Dict, save_dir::String="cluster_correction_results_complex")
    """
    Save results to a JLD2 file in the specified directory.
    """
    if !isdir(save_dir)
        mkpath(save_dir)
        println("📁 Created directory: $save_dir")
    end
    
    N, w, η = results["N"], results["w"], results["eta"]
    timestamp = replace(string(now()), ":" => "-")  # Make filename safe
    
    filename = "cluster_results_N$(N)_w$(w)_eta$(round(η, digits=3))_$(timestamp).jld2"
    filepath = joinpath(save_dir, filename)
    
    save(filepath, "results", results)
    println("💾 Results saved to: $filepath")
    
    return filepath
end

# Configuration
w_list = [4, 6, 8, 10]  # Available cluster weights
N_list = [5, 10]    # Available system sizes  
η_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # η values
nsamples = 1      # Number of samples per combination

# Get task ID from command line argument or environment variable
task_id = nothing
if length(ARGS) > 0
    task_id = parse(Int, ARGS[1])
elseif haskey(ENV, "SLURM_ARRAY_TASK_ID")
    task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
elseif haskey(ENV, "PBS_ARRAYID") 
    task_id = parse(Int, ENV["PBS_ARRAYID"])
else
    println("❌ No task ID provided!")
    println("Usage: julia cluster_exp_run.jl <task_id>")
    println("   or set SLURM_ARRAY_TASK_ID or PBS_ARRAYID environment variable")
    
    total_tasks = length(N_list) * length(w_list) * length(η_list)
    println("\nParameter combinations:")
    println("  N values: $N_list")
    println("  w values: $w_list")  
    println("  η values: $η_list")
    println("  Total combinations: $total_tasks")
    println("  Task IDs should range from 1 to $total_tasks")
    
    exit(1)
end

# Get parameters for this task
N, w, η = get_parameter_combination(task_id, N_list, w_list, η_list)

println("🎯 Task ID: $task_id")
println("📋 Parameters: N=$N, w=$w, η=$η")
println("🔢 Running $nsamples samples")

# Run the analysis
results = run_cluster_correction_analysis(N, w, η, nsamples)

if results !== nothing
    # Save results
    filepath = save_results(results)
    println("🎉 Task $task_id completed successfully!")
else
    println("❌ Task $task_id failed!")
    exit(1)
end