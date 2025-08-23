include("randompeps.jl")
using JLD2, FileIO, Dates, Statistics


function run_loop_correction_analysis(N::Int, w::Int, η::Float64, nsamples::Int)
    """
    Run naive loop correction analysis for given parameters.
    Computes free energy per site errors: log(Z) / (N*T).
    Returns Dict with results including average free energy per site errors.
    """
    
    println("="^80)
    println("🔄 Running Naive Loop Correction Analysis")
    println("Parameters: N=$N, w=$w, η=$η, nsamples=$nsamples")
    println("="^80)
    
    # Load loop data once
    cluster_data = nothing
    all_loops = nothing
    
    try
        cluster_data = load_latest_cluster_file(N, w)
        loop_objects = cluster_data.all_loops  
        all_loops = [loop_object.edges for loop_object in loop_objects]
        println("✅ Loaded loop data successfully ($(length(all_loops)) loops)")
    catch e
        println("❌ Failed to load loop data for N=$N, w=$w: $e")
        return nothing
    end
    
    # Storage for results
    bp_errors = Float64[]           # Free energy per site errors
    loop_errors = Float64[]         # Free energy per site errors with loop correction
    bp_fe_per_sites = Float64[]     # BP free energy per site values
    exact_fe_per_sites = Float64[]  # Exact free energy per site values
    loop_corrections = Float64[]    # Loop corrections (real parts)
    
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
            exact_PF = contract_peps_no_phys(peps; cutoff=1E-12, maxdim=2^N)
            exact_FE_per_site = log(Complex(exact_PF)) / (N*T)  # Complex log for safety
            
            # BP calculation
            adj_mat, edges, links = BP.get_adj_mat(tensors)
            messages = BP.get_messages(tensors, edges, links; random_part=0.01)
            messages = BP.message_passing(tensors, messages, edges, adj_mat; α=annealing, max_iters=maxiter)
            
            Z_l = BP.get_fixed_point_list(tensors, messages, adj_mat)
            
            # Compute BP free energy per site using complex arithmetic
            Z_l_complex = Complex.(Z_l)
            bp_FE_per_site = sum(log.(Z_l_complex)) / (N*T)
            T_normalized = BP.normalize_tensors(tensors, Z_l)
            
            # Naive loop correction - compute correction term to free energy per site
            loop_correction_sum = 0.0 + 0.0im  # Start with complex number for robustness
            valid_loops = 0
            total_loops = 0
            
            for loop in all_loops
                total_loops += 1
                if all([e in edges for e in loop])
                    try
                        loop_contri = scalar(BP.loop_contribution(loop, messages, T_normalized, edges, links, adj_mat))
                        # Add loop contribution to correction sum
                        loop_correction_sum += Complex(loop_contri)
                        valid_loops += 1
                    catch e
                        # Skip loops that cause errors
                        continue
                    end
                end
            end
            
            # Apply loop correction to free energy per site: log(Z + correction) / (N*T)
            loop_FE_per_site_correction = log(1.0 + loop_correction_sum)
            corrected_FE_per_site = bp_FE_per_site + loop_FE_per_site_correction
            
            # Debug output for first few samples
            if sample <= 3
                println("  Sample $sample debug:")
                println("    Total loops: $total_loops")
                println("    Valid loops: $valid_loops")
                println("    Exact FE per site: $exact_FE_per_site")
                println("    BP FE per site: $bp_FE_per_site")
                println("    Loop correction sum: $loop_correction_sum")
                println("    FE per site correction: $loop_FE_per_site_correction")
                println("    Corrected FE per site: $corrected_FE_per_site")
                println("    BP FE error: $(abs(bp_FE_per_site - exact_FE_per_site))")
                println("    Loop corrected FE error: $(abs(corrected_FE_per_site - exact_FE_per_site))")
            end
            
            # Calculate free energy per site errors
            bp_error = abs(bp_FE_per_site - exact_FE_per_site)
            loop_error = abs(corrected_FE_per_site - exact_FE_per_site)
            
            # Store results
            push!(bp_errors, bp_error)
            push!(loop_errors, loop_error)
            push!(bp_fe_per_sites, real(bp_FE_per_site))
            push!(exact_fe_per_sites, real(exact_FE_per_site))
            push!(loop_corrections, real(loop_FE_per_site_correction))  # Store FE correction
            
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
        "loop_fe_error_mean" => mean(loop_errors),
        "loop_fe_error_std" => std(loop_errors),
        "fe_improvement_mean" => mean(bp_errors) - mean(loop_errors),
        "fe_improvement_std" => std(bp_errors .- loop_errors),
        "loop_correction_mean" => mean(loop_corrections),
        "loop_correction_std" => std(loop_corrections),
        "timestamp" => string(now()),
        # Store all raw data for detailed analysis
        "bp_fe_errors" => bp_errors,
        "loop_fe_errors" => loop_errors,
        "bp_fe_per_sites" => bp_fe_per_sites,
        "exact_fe_per_sites" => exact_fe_per_sites,
        # Keep old names for backward compatibility
        "bp_errors" => bp_errors,
        "loop_errors" => loop_errors,
        "exact_fe_per_sites" => exact_fe_per_sites,
        "loop_corrections" => loop_corrections
    )
    
    println("✅ Analysis complete!")
    println("  Successful samples: $(length(bp_errors))/$nsamples")
    println("  BP free energy per site error (mean ± std): $(mean(bp_errors)) ± $(std(bp_errors))")
    println("  Loop free energy per site error (mean ± std): $(mean(loop_errors)) ± $(std(loop_errors))")
    println("  Free energy per site improvement: $(mean(bp_errors) - mean(loop_errors))")
    
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

function save_results(results::Dict, save_dir::String="loop_correction_results")
    """
    Save results to a JLD2 file in the specified directory.
    """
    if !isdir(save_dir)
        mkpath(save_dir)
        println("📁 Created directory: $save_dir")
    end
    
    N, w, η = results["N"], results["w"], results["eta"]
    timestamp = replace(string(now()), ":" => "-")  # Make filename safe
    
    filename = "loop_results_N$(N)_w$(w)_eta$(round(η, digits=3))_$(timestamp).jld2"
    filepath = joinpath(save_dir, filename)
    
    save(filepath, "results", results)
    println("💾 Results saved to: $filepath")
    
    return filepath
end

# Configuration
w_list = [4, 6, 8, 10]  # Available loop weights
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
    println("Usage: julia loop_corr_run.jl <task_id>")
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
results = run_loop_correction_analysis(N, w, η, nsamples)

if results !== nothing
    # Save results
    filepath = save_results(results)
    println("🎉 Task $task_id completed successfully!")
else
    println("❌ Task $task_id failed!")
    exit(1)
end
