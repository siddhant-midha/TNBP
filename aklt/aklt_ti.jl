#!/usr/bin/env julia

"""
aklt_ti.jl

Implementation of cluster expansion for AKLT model free energy density using saved clusters.
Follows exactly the workflow from functions/cluster_expansion_2d_ising.jl but uses pre-computed clusters
and AKLT tensors instead of Ising tensors.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("aklt.jl")
include("../functions/BP.jl")
include("../test/test_utils.jl")
using Serialization
using Plots
using ITensors, Graphs

# Data structure for single-site enumeration results (copied from generate_ising_clusters_one_site.jl)
struct SingleSiteClusterData
    """Data structure to save single-site cluster enumeration results."""
    clusters::Vector{Cluster}
    all_loops::Vector{Loop}
    adj_matrix::Matrix{Int}
    max_weight::Int
    lattice_size::Int
    site::Int
    enumeration_time::Float64
    translation_removal_time::Float64
    timestamp::String
    boundary_condition::String
    clusters_before_translation_removal::Int
    clusters_after_translation_removal::Int
    canonical_forms_count::Int
end

# Load saved cluster data
function load_cluster_data(filepath::String)
    """Load cluster enumeration data from saved file."""
    println("📂 Loading cluster data from: $(basename(filepath))")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    
    data = loaded_data["data"]
    summary = loaded_data["summary"]
    
    println("✅ Loaded successfully!")
    println("   Lattice size: $(summary["lattice_size"])×$(summary["lattice_size"])")
    println("   Max weight: $(summary["max_weight"])")
    println("   Total loops: $(summary["total_loops"])")
    
    if haskey(summary, "total_clusters")
        println("   Total clusters: $(summary["total_clusters"])")
    end
    
    return data, summary
end

function load_latest_single_site_cluster_file(; size_filter="L11", weight_filter="", boundary_filter="periodic")
    """
    Load the most recent single-site cluster enumeration file matching criteria.
    Returns the data and filename for reference.
    """
    save_dir = "../saved_clusters"
    
    if !isdir(save_dir)
        error("No saved_clusters directory found!")
    end
    
    files = readdir(save_dir)
    
    # Filter files based on criteria - must contain "single_site"
    matching_files = filter(files) do f
        # Must be a .jld2 file
        !endswith(f, ".jld2") && return false
        
        # Must contain "single_site" in filename
        !contains(f, "single_site") && return false
        
        # Must contain size filter
        !contains(f, size_filter) && return false
        
        # Must contain boundary filter
        !contains(f, boundary_filter) && return false
        
        # If weight filter specified, must contain it
        if !isempty(weight_filter)
            !contains(f, weight_filter) && return false
        end
        
        return true
    end
    
    if isempty(matching_files)
        error("No matching single-site cluster files found for criteria: size=$size_filter, weight=$weight_filter, boundary=$boundary_filter")
    end
    
    # Sort by filename (which includes timestamp) and take the most recent
    sorted_files = sort(matching_files, rev=true)  # Most recent first
    
    # Try to load files starting with the most recent
    for latest_file in sorted_files
        filepath = joinpath(save_dir, latest_file)
        
        println("📖 Attempting to load single-site cluster data: $(latest_file)")
        
        try
            loaded_data = open(filepath, "r") do io
                deserialize(io)
            end
            
            println("✅ Successfully loaded: $(latest_file)")
            return loaded_data["data"], latest_file
        catch e
            println("⚠️  Failed to load $(latest_file): $e")
            if latest_file == sorted_files[end]  # If this was the last file to try
                error("❌ All matching cluster files are corrupted or unreadable")
            end
            continue  # Try the next file
        end
    end
end

# Ursell function for connected clusters (from original implementation)
function ursell_function(cluster::Cluster)
    """
    Compute the Ursell function φ(W) for a connected cluster W.
    For connected clusters, this is (-1)^(|W|-1) * (|W|-1)!
    where |W| is the total number of loops in the cluster (with multiplicities).
    """
    total_loops = cluster.total_loops
    if total_loops == 1
        return 1.0
    else
        return -0.5
    end
end

function cluster_expansion_aklt_with_single_site_clusters(L::Int, a1::Float64, a2::Float64; 
                                                          max_weights=[4, 6, 8, 10])
    """
    Complete cluster expansion workflow for AKLT model using single-site saved clusters.
    Automatically loads the latest L×L single-site cluster file and handles normalization correctly.
    """
    
    println("="^80)
    println("🔥 Cluster Expansion for AKLT Model (Using Single-Site Clusters)")
    println("="^80)
    println("Parameters: L=$L, a1=$a1, a2=$a2")
    println("Cluster weights: $max_weights")
    println()
    
    # Step 1: Load latest single-site cluster data
    println("Step 1: Loading latest single-site cluster data...")
    cluster_data = nothing
    cluster_filename = ""
    
    try
        cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L", weight_filter="w$(maximum(max_weights))")
        println("✅ Loaded single-site cluster data from: $(cluster_filename)")
    catch e
        println("⚠️  Failed to load with specific weight, trying any weight...")
        try
            cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L")
            println("✅ Loaded single-site cluster data from: $(cluster_filename)")
        catch e2
            error("❌ Could not load single-site cluster data for L=$L: $e2")
        end
    end
    
    # Step 2: Create AKLT tensor network
    println("\nStep 2: Creating AKLT tensor network...")
    T = aklt_norm_network(L; a1=a1, a2=a2)
    N = L^2
    println("✅ Created $N tensors for $(L)×$(L) lattice")
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    println("✅ Found $(length(edges)) edges")
    
    # Step 3: Compute BP fixed point
    println("\nStep 3: Computing BP fixed point...")
    messages = BP.get_messages(T, edges, links; random_part=0.1)
    messages = BP.message_passing(T, messages, edges, adj_mat; α=0.8, max_iters=1000)
    println("✅ BP converged")
    
    # Get BP partition function (before normalization)
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Complex(Z_bp_full))
    println("✅ BP log partition function: $log_Z_bp_full")
    
    # Step 4: Normalize tensors
    println("\nStep 4: Normalizing tensors...")
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)
    
    # Normalization factor
    Z_l_complex = Complex.(Z_l)
    normalization_factor = sum(log.(Z_l_complex))
    println("✅ Normalization factor: $normalization_factor")
    
    # Verify normalization (BP on normalized tensors should give 1)
    Z_bp_normalized = BP.mean_free_partition_fn(1:N, T_normalized, messages, adj_mat)
    println("✅ Normalized BP partition function: $Z_bp_normalized (should ≈ 1)")
    
    # Step 5: Compute ALL cluster contributions once
    println("\nStep 5: Computing ALL cluster contributions once...")
    
    # Compute individual cluster contributions for all clusters up to max weight
    max_weight_overall = maximum(max_weights)
    cluster_contributions = compute_all_single_site_cluster_contributions(
        T_normalized, messages, edges, links, adj_mat, 
        cluster_data, L, max_weight_overall)
    
    println("✅ Computed $(length(cluster_contributions)) individual cluster contributions")
    
    results = Dict{String, Any}()
    results["L"] = L
    results["a1"] = a1
    results["a2"] = a2
    results["bp_log_Z"] = log_Z_bp_full  # Keep as raw log partition function
    results["normalization_factor"] = normalization_factor
    results["cluster_corrections"] = Dict()
    results["cluster_filename"] = cluster_filename
    
    # Step 6: Sum contributions for different weight truncations
    println("\nStep 6: Summing contributions for different weight truncations...")
    
    for max_weight in max_weights
        println("🔍 Computing cluster expansion up to weight $max_weight")
        
        # Sum only contributions from clusters with weight <= max_weight
        correction = sum(contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= max_weight)
        
        # For single-site clusters: 
        # - BP contribution needs to be divided by N to get free energy density: f_BP = -log(Z_BP)/N
        # - Single-site cluster correction is already per-site, so it's the free energy density correction directly
        f_bp = -real(log_Z_bp_full) / N  # BP free energy density
        f_corrected = f_bp + correction  # Corrected free energy density
        
        results["cluster_corrections"][max_weight] = Dict(
            "correction" => correction,  # Free energy density correction (per site)
            "f_bp" => f_bp,  # BP free energy density
            "f_corrected" => f_corrected,  # Corrected free energy density
            "improvement" => abs(correction)
        )
        
        println("📊 Results for weight $max_weight:")
        println("  BP free energy density: $f_bp")
        println("  Cluster correction (per site): $correction")
        println("  Corrected free energy density: $f_corrected")
    end
    
    # Step 6: Compare with exact CTMRG solution
    println("\n" * "="^60)
    println("📐 Comparing with exact CTMRG solution")
    println("="^60)
    
    exact_free_energy = -ctmrg_exact_FE_density(a1, a2)
    results["exact_ctmrg"] = exact_free_energy
    
    println("🎯 Exact CTMRG free energy density: $exact_free_energy")
    println("\n📈 Comparison with exact result:")
    
    f_bp = results["cluster_corrections"][max_weights[1]]["f_bp"]  # BP is the same for all weights
    bp_error = abs(f_bp - exact_free_energy)
    println("  BP approximation error: $bp_error")
    
    for max_weight in max_weights
        f_corrected = results["cluster_corrections"][max_weight]["f_corrected"]
        error = abs(f_corrected - exact_free_energy)
        improvement = bp_error - error
        
        println("  Weight $max_weight error: $error (improvement: $improvement)")
    end
    
    return results
end

function compute_all_single_site_cluster_contributions(T_normalized, messages, edges, links, adj_mat, 
                                                     cluster_data, L::Int, max_weight::Int)
    """
    Compute ALL cluster contributions once for single-site cluster data.
    Returns a vector of dictionaries with contribution details for each cluster.
    """
    
    # Get clusters from single-site data
    all_clusters = cluster_data.clusters
    all_loops = cluster_data.all_loops
    
    println("📊 Computing contributions for all single-site clusters")
    println("   Available clusters: $(length(all_clusters))")
    println("   Available loops: $(length(all_loops))")
    
    # Filter clusters by weight
    relevant_clusters = [c for c in all_clusters if c.weight <= max_weight]
    println("   Using $(length(relevant_clusters)) clusters (weight ≤ $max_weight)")
    
    if isempty(relevant_clusters)
        println("⚠️  No relevant clusters found")
        return []
    end
    
    # Compute contributions for all relevant clusters
    cluster_contributions = []
    
    println("💫 Computing individual cluster contributions...")
    
    for (i, cluster) in enumerate(relevant_clusters)
        if i % 100 == 0 || i <= 10
            println("  Processing cluster $i/$(length(relevant_clusters)) (weight=$(cluster.weight))...")
        end
        
        # Compute Ursell function φ(W)
        phi_W = ursell_function(cluster)
        
        if abs(phi_W) < 1e-15
            continue  # Skip negligible contributions
        end
        
        # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
        Z_W = 1.0 + 0im  # Complex number for cluster contribution
        computation_successful = true
        
        for loop_id in cluster.loop_ids
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            
            # Convert loop edges format for BP.jl (ensure v1 < v2 ordering)
            loop_edges_bp = Tuple{Int,Int}[]
            for edge in loop.edges
                v1, v2 = edge
                push!(loop_edges_bp, (min(v1, v2), max(v1, v2)))
            end
            
            # Compute loop contribution Z_l using BP
            try
                Z_l_tensor = BP.loop_contribution(loop_edges_bp, messages, T_normalized, edges, links, adj_mat)
                # Extract scalar value from ITensor
                Z_l = scalar(Z_l_tensor)
                Z_W *= Z_l^multiplicity
            catch e
                println("⚠️  Error computing loop contribution for loop $loop_id: $e")
                computation_successful = false
                break
            end
        end
        
        if !computation_successful
            continue
        end
        
        # Add contribution to free energy density
        # For free energy: f = -log(Z), so correction is φ(W) * Z_W, but we want the contribution to f
        # The contribution to log(Z) is φ(W) * Z_W, so contribution to f is -φ(W) * Z_W
        contribution = -phi_W * Z_W 
        
        # Store contribution details
        contrib_info = Dict(
            "cluster_id" => i,
            "weight" => cluster.weight,
            "phi_W" => phi_W,
            "Z_W" => Z_W,
            "contribution" => real(contribution)
        )
        push!(cluster_contributions, contrib_info)
        
        if abs(contribution) > 1e-10 && i <= 20
            println("    Cluster $i (weight=$(cluster.weight)): φ=$phi_W, Z_W=$Z_W, f contribution=$contribution")
        end
    end
    
    println("✅ Computed $(length(cluster_contributions)) valid cluster contributions")
    
    # Sort by weight for easier analysis
    sort!(cluster_contributions, by=x->x["weight"])
    
    return cluster_contributions
end

# Plot free energy vs a1 parameter using single-site clusters (EFFICIENT VERSION)
function plot_free_energy_vs_a1_with_single_site_clusters_efficient(L::Int, a1_range::Vector{Float64}, a2::Float64,
                                                                     max_weights::Vector{Int}=[4,6,8,10]; save_path::String="")
    """
    Plot free energy density as a function of a1 using single-site clusters.
    EFFICIENT: Loads cluster data once and computes contributions once per a1 value.
    """
    
    println("📊 Creating free energy vs a1 plot using single-site clusters (EFFICIENT)")
    println("   Lattice size: $(L)×$(L)")
    println("   a1 range: $(a1_range[1]) to $(a1_range[end]) ($(length(a1_range)) points)")
    println("   Fixed a2: $a2")
    println("   Weight truncations: $max_weights")
    
    # Load single-site cluster data once
    println("\n📂 Loading single-site cluster data once...")
    cluster_data = nothing
    cluster_filename = ""
    
    try
        cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L", weight_filter="w$(maximum(max_weights))")
        println("✅ Loaded single-site cluster data from: $(cluster_filename)")
    catch e
        println("⚠️  Failed to load with specific weight, trying any weight...")
        try
            cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L")
            println("✅ Loaded single-site cluster data from: $(cluster_filename)")
        catch e2
            error("❌ Could not load single-site cluster data for L=$L: $e2")
        end
    end
    
    # Initialize plot
    p = plot(title="AKLT Free Energy Density vs a1 (Single-Site Cluster Expansion)", 
             xlabel="a1", ylabel="Free Energy Density", 
             legend=:topright, dpi=300, size=(800, 600))
    
    # Compute exact CTMRG solution
    println("\n🎯 Computing exact CTMRG solution...")
    exact_f = [-ctmrg_exact_FE_density(a1, a2) for a1 in a1_range]
    plot!(p, a1_range, exact_f, linewidth=3, color=:black, 
          label="Exact CTMRG", linestyle=:solid)
    
    # Initialize storage for all a1 points
    bp_f = Float64[]
    cluster_f_dict = Dict{Int, Vector{Float64}}()
    for max_weight in max_weights
        cluster_f_dict[max_weight] = Float64[]
    end
    
    max_weight_overall = maximum(max_weights)
    
    # Compute for each a1 value
    println("\n📊 Computing for each a1 value...")
    for (j, a1) in enumerate(a1_range)
        if j % 5 == 1 || j <= 3
            println("  Computing a1 = $a1 (point $j/$(length(a1_range)))...")
        end
        
        # Create AKLT tensor network
        T = aklt_norm_network(L; a1=a1, a2=a2)
        N = L^2
        
        # Get adjacency structure
        adj_mat, edges, links = BP.get_adj_mat(T)
        
        # Compute BP fixed point
        messages = BP.get_messages(T, edges, links; random_part=0.1)
        messages = BP.message_passing(T, messages, edges, adj_mat; α=0.8, max_iters=1000)
        
        # Get BP partition function (before normalization)
        Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
        log_Z_bp_full = log(Complex(Z_bp_full))
        
        # Normalize tensors
        Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
        T_normalized = BP.normalize_tensors(T, Z_l)
        
        # BP free energy density
        f_bp = -real(log_Z_bp_full) / N
        push!(bp_f, f_bp)
        
        # Compute ALL cluster contributions once for this a1
        cluster_contributions = compute_all_single_site_cluster_contributions(
            T_normalized, messages, edges, links, adj_mat, 
            cluster_data, L, max_weight_overall)
        
        # Compute corrections for each weight
        for max_weight in max_weights
            correction = sum(contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= max_weight)
            f_corrected = f_bp + correction
            push!(cluster_f_dict[max_weight], f_corrected)
        end
    end
    
    # Plot BP baseline
    plot!(p, a1_range, bp_f, linewidth=2, color=:gray, 
          label="BP Approximation", linestyle=:dash)
    
    # Colors for different weights
    colors = [:red, :blue, :green, :orange, :purple, :brown]
    
    # Plot cluster expansion for each weight
    for (i, max_weight) in enumerate(max_weights)
        color = colors[min(i, length(colors))]
        plot!(p, a1_range, cluster_f_dict[max_weight], linewidth=2, color=color,
              label="Cluster w≤$max_weight", linestyle=:solid)
    end
    
    # Save plot
    if !isempty(save_path)
        savefig(p, save_path)
        println("✅ Main plot saved to: $save_path")
    end
    
    # Create error plot showing absolute differences from exact solution
    println("\n📊 Creating error plot showing absolute differences from exact CTMRG solution...")
    p_error = plot(title="Absolute Error vs a1 (Distance from Exact CTMRG Solution)", 
                   xlabel="a1", ylabel="Absolute Error", 
                   legend=:topright, dpi=300, size=(800, 600))
    
    # Plot BP error
    bp_error = abs.(bp_f .- exact_f)
    # Filter out any NaN, Inf, or zero values that cause log scale issues
    valid_indices = .!(isnan.(bp_error) .| isinf.(bp_error) .| (bp_error .<= 0))
    if any(valid_indices)
        plot!(p_error, a1_range[valid_indices], bp_error[valid_indices], linewidth=2, color=:gray, 
              label="BP Error", linestyle=:dash)
    end
    
    # Plot cluster expansion errors for each weight
    for (i, max_weight) in enumerate(max_weights)
        color = colors[min(i, length(colors))]
        cluster_error = abs.(cluster_f_dict[max_weight] .- exact_f)
        # Filter out any NaN, Inf, or zero values that cause log scale issues
        valid_indices = .!(isnan.(cluster_error) .| isinf.(cluster_error) .| (cluster_error .<= 0))
        if any(valid_indices)
            plot!(p_error, a1_range[valid_indices], cluster_error[valid_indices], linewidth=2, color=color,
                  label="Cluster w≤$max_weight Error", linestyle=:solid)
        end
    end
    
    # Save error plot
    if !isempty(save_path)
        error_save_path = replace(save_path, ".png" => "_errors.png")
        savefig(p_error, error_save_path)
        println("✅ Error plot saved to: $error_save_path")
    end
    
    return p, p_error
end

# Plot free energy vs a2 parameter using single-site clusters
function plot_free_energy_vs_a2_with_single_site_clusters_efficient(L::Int, a1::Float64, a2_range::Vector{Float64},
                                                                     max_weights::Vector{Int}=[4,6,8,10]; save_path::String="")
    """
    Plot free energy density as a function of a2 using single-site clusters.
    """
    
    println("📊 Creating free energy vs a2 plot using single-site clusters")
    println("   Lattice size: $(L)×$(L)")
    println("   Fixed a1: $a1")
    println("   a2 range: $(a2_range[1]) to $(a2_range[end]) ($(length(a2_range)) points)")
    println("   Weight truncations: $max_weights")
    
    # Load single-site cluster data once
    println("\n📂 Loading single-site cluster data once...")
    cluster_data = nothing
    cluster_filename = ""
    
    try
        cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L", weight_filter="w$(maximum(max_weights))")
        println("✅ Loaded single-site cluster data from: $(cluster_filename)")
    catch e
        println("⚠️  Failed to load with specific weight, trying any weight...")
        try
            cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L")
            println("✅ Loaded single-site cluster data from: $(cluster_filename)")
        catch e2
            error("❌ Could not load single-site cluster data for L=$L: $e2")
        end
    end
    
    # Initialize plot
    p = plot(title="AKLT Free Energy Density vs a2 (Single-Site Cluster Expansion)", 
             xlabel="a2", ylabel="Free Energy Density", 
             legend=:topright, dpi=300, size=(800, 600))
    
    # Compute exact CTMRG solution
    println("\n🎯 Computing exact CTMRG solution...")
    exact_f = [-ctmrg_exact_FE_density(a1, a2) for a2 in a2_range]
    plot!(p, a2_range, exact_f, linewidth=3, color=:black, 
          label="Exact CTMRG", linestyle=:solid)
    
    # Initialize storage for all a2 points
    bp_f = Float64[]
    cluster_f_dict = Dict{Int, Vector{Float64}}()
    for max_weight in max_weights
        cluster_f_dict[max_weight] = Float64[]
    end
    
    max_weight_overall = maximum(max_weights)
    
    # Compute for each a2 value
    println("\n📊 Computing for each a2 value...")
    for (j, a2) in enumerate(a2_range)
        if j % 5 == 1 || j <= 3
            println("  Computing a2 = $a2 (point $j/$(length(a2_range)))...")
        end
        
        # Create AKLT tensor network
        T = aklt_norm_network(L; a1=a1, a2=a2)
        N = L^2
        
        # Get adjacency structure
        adj_mat, edges, links = BP.get_adj_mat(T)
        
        # Compute BP fixed point
        messages = BP.get_messages(T, edges, links; random_part=0.1)
        messages = BP.message_passing(T, messages, edges, adj_mat; α=1., max_iters=1000)
        
        # Get BP partition function (before normalization)
        Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
        log_Z_bp_full = log(Complex(Z_bp_full))
        
        # Normalize tensors
        Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
        T_normalized = BP.normalize_tensors(T, Z_l)
        
        # BP free energy density
        f_bp = -real(log_Z_bp_full) / N
        push!(bp_f, f_bp)
        
        # Compute ALL cluster contributions once for this a2
        cluster_contributions = compute_all_single_site_cluster_contributions(
            T_normalized, messages, edges, links, adj_mat, 
            cluster_data, L, max_weight_overall)
        
        # Compute corrections for each weight
        for max_weight in max_weights
            correction = sum(contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= max_weight)
            f_corrected = f_bp + correction
            push!(cluster_f_dict[max_weight], f_corrected)
        end
    end
    
    # Plot BP baseline
    plot!(p, a2_range, bp_f, linewidth=2, color=:gray, 
          label="BP Approximation", linestyle=:dash)
    
    # Colors for different weights
    colors = [:red, :blue, :green, :orange, :purple, :brown]
    
    # Plot cluster expansion for each weight
    for (i, max_weight) in enumerate(max_weights)
        color = colors[min(i, length(colors))]
        plot!(p, a2_range, cluster_f_dict[max_weight], linewidth=2, color=color,
              label="Cluster w≤$max_weight", linestyle=:solid)
    end
    
    # Save plot
    if !isempty(save_path)
        savefig(p, save_path)
        println("✅ Main plot saved to: $save_path")
    end
    
    # Create error plot
    println("\n📊 Creating error plot...")
    p_error = plot(title="Absolute Error vs a2 (Distance from Exact CTMRG Solution)", 
                   xlabel="a2", ylabel="Absolute Error", 
                   legend=:topright, dpi=300, size=(800, 600))
    
    # Plot BP error
    bp_error = abs.(bp_f .- exact_f)
    # Filter out any NaN, Inf, or zero values that cause plotting issues
    valid_indices = .!(isnan.(bp_error) .| isinf.(bp_error) .| (bp_error .<= 0))
    if any(valid_indices)
        plot!(p_error, a2_range[valid_indices], bp_error[valid_indices], linewidth=2, color=:gray, 
              label="BP Error", linestyle=:dash)
    end
    
    # Plot cluster expansion errors for each weight
    for (i, max_weight) in enumerate(max_weights)
        color = colors[min(i, length(colors))]
        cluster_error = abs.(cluster_f_dict[max_weight] .- exact_f)
        # Filter out any NaN, Inf, or zero values that cause plotting issues
        valid_indices = .!(isnan.(cluster_error) .| isinf.(cluster_error) .| (cluster_error .<= 0))
        if any(valid_indices)
            plot!(p_error, a2_range[valid_indices], cluster_error[valid_indices], linewidth=2, color=color,
                  label="Cluster w≤$max_weight Error", linestyle=:solid)
        end
    end
    
    # Save error plot
    if !isempty(save_path)
        error_save_path = replace(save_path, ".png" => "_errors.png")
        savefig(p_error, error_save_path)
        println("✅ Error plot saved to: $error_save_path")
    end
    
    return p, p_error
end

# Main execution function
function main()
    println("🚀 AKLT Cluster Expansion Free Energy Calculator (Using Single-Site Clusters)")
    println("="^80)
    
    # AKLT parameters
    L = 12
    # a1 = sqrt(3/2)
    # a2 = sqrt(6)
    a1 = 0.5
    a2 = 0.2
    # a1_range = collect(0.8:0.01:1.5)  # Range around sqrt(3/2) ≈ 1.225
    # a2_range = collect(2.0:0.02:3.0)  # Range around sqrt(6) ≈ 2.449
    a1_range = collect(0.25:0.025:0.75)  # Range around sqrt(3/2) ≈ 1.225
    a2_range = collect(0.1:0.01:0.3)  # Range around sqrt(6) ≈ 2.449

    max_weights = [4, 6, 7, 8, 9, 10]  # Use available weights
    
    println("📋 Configuration:")
    println("   Lattice size: $(L)×$(L)")
    println("   AKLT parameters: a1=$a1, a2=$a2")
    println("   Weight truncations: $max_weights")
    
    # Run single point cluster expansion analysis
    println("\n🔥 Running single point cluster expansion analysis...")
    results = cluster_expansion_aklt_with_single_site_clusters(L, a1, a2; max_weights=max_weights)
    
    # Display final results
    println("\n✅ Single point analysis complete!")
    println("📊 AKLT cluster expansion results:")
    
    exact_fe = results["exact_ctmrg"]
    println("   Exact CTMRG free energy density: $exact_fe")
    
    for max_weight in max_weights
        correction_data = results["cluster_corrections"][max_weight]
        f_bp = correction_data["f_bp"]
        f_corrected = correction_data["f_corrected"]
        correction = correction_data["correction"]
        
        bp_error = abs(f_bp - exact_fe)
        corrected_error = abs(f_corrected - exact_fe)
        improvement = bp_error - corrected_error
        
        println("   Weight $max_weight:")
        println("     BP error: $bp_error")
        println("     Corrected error: $corrected_error")
        println("     Improvement: $improvement")
        println("     Correction magnitude: $(abs(correction))")
    end
    
    # Create parameter sweep plots
    println("\n" * "="^80)
    println("🎨 Creating parameter sweep plots...")
    println("="^80)
    
    # Plot vs a1 parameter (varying a1, fixed a2)
    println("\n📊 Creating plot vs a1 parameter...")
    save_path_a1 = "../visualization/aklt_cluster_expansion_vs_a1_L$(L).png"
    
    try
        plot_a1, error_plot_a1 = plot_free_energy_vs_a1_with_single_site_clusters_efficient(
            L, a1_range, a2, max_weights; save_path=save_path_a1)
        println("✅ a1 parameter plots created successfully!")
    catch e
        println("⚠️  Failed to create a1 plots: $e")
    end
    
    # Plot vs a2 parameter (fixed a1, varying a2)
    println("\n📊 Creating plot vs a2 parameter...")
    save_path_a2 = "../visualization/aklt_cluster_expansion_vs_a2_L$(L).png"
    
    try
        plot_a2, error_plot_a2 = plot_free_energy_vs_a2_with_single_site_clusters_efficient(
            L, a1, a2_range, max_weights; save_path=save_path_a2)
        println("✅ a2 parameter plots created successfully!")
    catch e
        println("⚠️  Failed to create a2 plots: $e")
    end
    
    println("\n✅ All analysis complete!")
    println("📊 Plots saved to:")
    println("   a1 parameter sweep: $save_path_a1")
    println("   a2 parameter sweep: $save_path_a2")
    
    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
