include("randompeps.jl")
include("cluster_exp.jl")
using JLD2, FileIO, Dates, Statistics


function run_cluster_correction_analysis(N::Int, w::Int, η::Float64, nsamples::Int)
    """
    Run cluster correction analysis for given parameters.
    Returns Dict with results including average errors.
    """
    
    println("="^80)
    println("🔥 Running Cluster Correction Analysis")
    println("Parameters: N=$N, w=$w, η=$η, nsamples=$nsamples")
    println("="^80)
    
    # Load cluster data once
    cluster_data = nothing
    clusters_by_site = nothing
    all_loops = nothing
    
    try
        cluster_data = load_latest_cluster_file(N, w)
        clusters_by_site = cluster_data.clusters_by_site
        all_loops = cluster_data.all_loops
        println("✅ Loaded cluster data successfully")
        println("🔍 Debug: clusters_by_site has $(length(clusters_by_site)) sites")
        println("🔍 Debug: all_loops has $(length(all_loops)) loops")
        if length(clusters_by_site) > 0
            first_site_clusters = length(clusters_by_site[1])
            println("🔍 Debug: First site has $first_site_clusters clusters")
        end
    catch e
        println("❌ Failed to load cluster data for N=$N, w=$w: $e")
        return nothing
    end
    
    # Storage for results
    bp_errors = Float64[]
    cluster_errors = Float64[]
    bp_fe_densities = Float64[]
    exact_fe_densities = Float64[]
    cluster_corrections = Float64[]
    
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
            tensors, peps = peps_controllable(N, T; η=η, ti=ti, orthog=orthog)
            exact_FE_density = (log(real(contract_peps_no_phys(peps; cutoff=1E-12, maxdim=2^N)))) / (N*T)
            
            # BP calculation
            adj_mat, edges, links = BP.get_adj_mat(tensors)
            messages = BP.get_messages(tensors, edges, links; random_part=0.0)
            messages = BP.message_passing(tensors, messages, edges, adj_mat; α=annealing, max_iters=maxiter)
            
            Z_l = BP.get_fixed_point_list(tensors, messages, adj_mat)
            bp_FE_density = sum(log.(real.(Z_l))) / (N*T)
            T_normalized = BP.normalize_tensors(tensors, Z_l)
            
            # Cluster correction (using manual deduplication)
            clustercorrx = 0.0
            
            # Deduplicate clusters manually (we know Set doesn't work)
            clusters_by_signature = Dict{Tuple, Vector{Tuple{Int, Cluster}}}()
            
            for site in 1:(N*T)
                clusters = clusters_by_site[site]
                for cluster in clusters
                    signature = (cluster.weight, cluster.total_loops, sort(cluster.loop_ids), sort(collect(cluster.multiplicities)))
                    if !haskey(clusters_by_signature, signature)
                        clusters_by_signature[signature] = []
                    end
                    push!(clusters_by_signature[signature], (site, cluster))
                end
            end
            
            # Compute contributions for unique clusters
            valid_contributions = 0
            total_contributions = 0
            for (sig, cluster_list) in clusters_by_signature
                cluster = cluster_list[1][2]  # Take first representative
                contribution = cluster_contr(T_normalized, messages, edges, links, adj_mat, cluster, all_loops)
                total_contributions += 1
                if !isnan(contribution) && isfinite(contribution)
                    clustercorrx += contribution
                    valid_contributions += 1
                end
            end
            
            cluster_FE_density_correction = clustercorrx / (N*T)
            
            # Debug output for first sample only (to avoid login node overload)
            if sample == 1
                println("  🔍 Sample 1 debug:")
                println("    Unique cluster signatures: $(length(clusters_by_signature))")
                println("    Total contributions attempted: $total_contributions")
                println("    Valid contributions: $valid_contributions")
                println("    Raw cluster correction sum: $clustercorrx")
                println("    Cluster FE correction density: $cluster_FE_density_correction")
                println("    BP error: $(abs(bp_FE_density - exact_FE_density))")
                println("    Will be cluster error: $(abs(bp_FE_density + cluster_FE_density_correction - exact_FE_density))")
            end
            
            # Calculate errors
            bp_error = abs(bp_FE_density - exact_FE_density)
            cluster_error = abs(bp_FE_density + cluster_FE_density_correction - exact_FE_density)
            
            # Store results
            push!(bp_errors, bp_error)
            push!(cluster_errors, cluster_error)
            push!(bp_fe_densities, bp_FE_density)
            push!(exact_fe_densities, exact_FE_density)
            push!(cluster_corrections, cluster_FE_density_correction)
            
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
        "bp_error_mean" => mean(bp_errors),
        "bp_error_std" => std(bp_errors),
        "cluster_error_mean" => mean(cluster_errors),
        "cluster_error_std" => std(cluster_errors),
        "improvement_mean" => mean(bp_errors) - mean(cluster_errors),
        "improvement_std" => std(bp_errors .- cluster_errors),
        "cluster_correction_mean" => mean(cluster_corrections),
        "cluster_correction_std" => std(cluster_corrections),
        "timestamp" => string(now()),
        # Store all raw data for detailed analysis
        "bp_errors" => bp_errors,
        "cluster_errors" => cluster_errors,
        "bp_fe_densities" => bp_fe_densities,
        "exact_fe_densities" => exact_fe_densities,
        "cluster_corrections" => cluster_corrections
    )
    
    println("✅ Analysis complete!")
    println("  Successful samples: $(length(bp_errors))/$nsamples")
    println("  BP error (mean ± std): $(mean(bp_errors)) ± $(std(bp_errors))")
    println("  Cluster error (mean ± std): $(mean(cluster_errors)) ± $(std(cluster_errors))")
    println("  Improvement: $(mean(bp_errors) - mean(cluster_errors))")
    
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

function save_results(results::Dict, save_dir::String="cluster_correction_results")
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
nsamples = 1000      # Number of samples per combination

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