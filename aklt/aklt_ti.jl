#!/usr/bin/env julia

"""
aklt_ti.jl

Implementation of cluster expansion for AKLT model free energy density using saved clusters.
Follows exactly the workflow from functions/cluster_expansion_2d_ising.jl but uses pre-computed clusters
and AKLT tensors instead of Ising tensors.
"""

include("aklt.jl")
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("../functions/BP.jl")
using Serialization
using Plots

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
    p_error = plot(legend=:bottomright, dpi=300, size=(800, 600), yscale=:log10)
    
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
        error_save_path = replace(save_path, ".pdf" => "_errors.pdf")
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
    p_error = plot(legend=:topright, dpi=300, size=(800, 600), yscale=:log10)
    
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
        error_save_path = replace(save_path, ".pdf" => "_errors.pdf")
        savefig(p_error, error_save_path)
        println("✅ Error plot saved to: $error_save_path")
    end
    
    return p, p_error
end

# Plot error vs weight
function plot_error_vs_weight(L::Int, a1::Float64, a2::Float64, max_weights::Vector{Int}; 
                               save_path::String="", fontfamily="Times New Roman")
    """
    Plot the absolute errors of cluster expansion as a function of cluster weight.
    
    This function computes the free energy using BP and cluster expansion for different
    maximum cluster weights, and plots the absolute errors relative to the exact CTMRG solution.
    
    Args:
        L: Linear size of the lattice
        a1, a2: AKLT model parameters
        max_weights: Vector of maximum weights to consider (should include weights in ascending order)
        save_path: Optional path to save the plot
        fontfamily: Font family to use for the plot
        
    Returns:
        The plot object
    """
    println("📊 Creating error vs weight plot for AKLT model")
    println("   Parameters: L=$L, a1=$a1, a2=$a2")
    println("   Weight truncations: $max_weights")
    
    # Load cluster data
    println("\n📂 Loading single-site cluster data...")
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
    
    # Create AKLT tensor network
    println("\n📊 Creating AKLT tensor network...")
    T = aklt_norm_network(L; a1=a1, a2=a2)
    N = L^2
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    
    # Compute BP fixed point
    println("📊 Computing BP fixed point...")
    messages = BP.get_messages(T, edges, links; random_part=0.1)
    messages = BP.message_passing(T, messages, edges, adj_mat; α=0.8, max_iters=1000)
    
    # Get BP partition function and free energy density
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Complex(Z_bp_full))
    f_bp = -real(log_Z_bp_full) / N
    
    # Normalize tensors
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)
    
    # Compute exact CTMRG solution
    exact_fe = -ctmrg_exact_FE_density(a1, a2)
    
    # Compute ALL cluster contributions once
    println("📊 Computing cluster contributions...")
    cluster_contributions = compute_all_single_site_cluster_contributions(
        T_normalized, messages, edges, links, adj_mat, 
        cluster_data, L, maximum(max_weights))
    
    # Compute errors for each weight
    weights = [0; max_weights]  # Include BP (weight 0)
    errors = zeros(length(weights))
    
    # BP error (weight 0)
    errors[1] = abs(f_bp - exact_fe)
    
    # Cluster expansion errors
    for (i, w) in enumerate(max_weights)
        correction = sum(contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= w)
        f_corrected = f_bp + correction
        errors[i+1] = abs(f_corrected - exact_fe)
    end
    
    # Create plot
    println("📊 Creating plot...")
    p = plot(legend=:topright, dpi=300, size=(800, 600),
             yscale=:log10, fontfamily=fontfamily, xticks=(weights, string.(weights)))
    
    # Plot errors
    plot!(p, weights, errors, linewidth=2, color=:blue, 
          marker=:circle, markersize=6, label="Error")
        
    # # Annotate the plot
    # annotate!(p, [(weights[end], errors[1], 
    #           text("BP", fontfamily=fontfamily, pointsize=10, :right, :top))])
    # annotate!(p, [(weights[end], errors[end], 
    #           text("w≤$(max_weights[end])", fontfamily=fontfamily, pointsize=10, :right, :bottom))])
    
    # Save plot
    if !isempty(save_path)
        savefig(p, save_path)
        println("✅ Plot saved to: $save_path")
    end
    
    return p
end

# Main execution function
function main()
    println("🚀 AKLT Cluster Expansion Free Energy Calculator (Using Single-Site Clusters)")
    println("="^80)
    
    # AKLT parameters
    L = 10
    # a1 = sqrt(3/2)
    # a2 = sqrt(6)
    a1 = 0.5
    a2 = 0.2
    # a1_range = collect(0.8:0.01:1.5)  # Range around sqrt(3/2) ≈ 1.225
    # a2_range = collect(2.0:0.02:3.0)  # Range around sqrt(6) ≈ 2.449
    a1_range = collect(0.25:0.025:0.75)  # Range around sqrt(3/2) ≈ 1.225
    a2_range = collect(0.1:0.01:0.3)  # Range around sqrt(6) ≈ 2.449

    max_weights = [4, 6, 7, 8, 9, 10]  # Use available weights
    results = cluster_expansion_aklt_with_single_site_clusters(L, a1, a2; max_weights=max_weights)
    
    
    exact_fe = results["exact_ctmrg"]
    
    for max_weight in max_weights
        correction_data = results["cluster_corrections"][max_weight]
        f_bp = correction_data["f_bp"]
        f_corrected = correction_data["f_corrected"]
        correction = correction_data["correction"]
        
        bp_error = abs(f_bp - exact_fe)
        corrected_error = abs(f_corrected - exact_fe)
        improvement = bp_error - corrected_error
    end
    
    # Create parameter sweep plots
    println("\n" * "="^80)
    println("🎨 Creating parameter sweep plots...")
    println("="^80)
    
    # Plot vs a1 parameter (varying a1, fixed a2)
    println("\n📊 Creating plot vs a1 parameter...")
    save_path_a1 = "../visualization/aklt_cluster_expansion_vs_a1_L$(L).pdf"
    
    try
        plot_a1, error_plot_a1 = plot_free_energy_vs_a1_with_single_site_clusters_efficient(
            L, a1_range, a2, max_weights; save_path=save_path_a1)
        println("✅ a1 parameter plots created successfully!")
    catch e
        println("⚠️  Failed to create a1 plots: $e")
    end
    
    # Plot vs a2 parameter (fixed a1, varying a2)
    println("\n📊 Creating plot vs a2 parameter...")
    save_path_a2 = "../visualization/aklt_cluster_expansion_vs_a2_L$(L).pdf"
    
    try
        plot_a2, error_plot_a2 = plot_free_energy_vs_a2_with_single_site_clusters_efficient(
            L, a1, a2_range, max_weights; save_path=save_path_a2)
        println("✅ a2 parameter plots created successfully!")
    catch e
        println("⚠️  Failed to create a2 plots: $e")
    end
    
    # Plot error vs weight
    println("\n📊 Creating error vs weight plot...")
    save_path_error = "../visualization/aklt_error_vs_weight_L$(L)_a1$(a1)_a2$(a2).pdf"
    
    try
        error_plot = plot_error_vs_weight(L, a1, a2, max_weights; 
                                         save_path=save_path_error,
                                         fontfamily="Times New Roman")
        println("✅ Error vs weight plot created successfully!")
    catch e
        println("⚠️  Failed to create error vs weight plot: $e")
    end
    
    println("\n✅ All analysis complete!")
    println("📊 Plots saved to:")
    println("   a1 parameter sweep: $save_path_a1")
    println("   a2 parameter sweep: $save_path_a2")
    println("   error vs weight: $save_path_error")
    
    return results
end


L = 20
max_weights = [4,6,8] 

# a1 = sqrt(3/2); a2 = sqrt(6)

a1 = 0.5
a2_range = vcat(collect(0.5:0.1:0.9), collect(0.9:0.0025:1.1), collect(1.2:0.1:1.9), collect(1.9:0.01:2.1),collect(2.2:0.1:3.0))#collect(0.5:0.01:3.0) 
save_path_a2 = "../visualization/aklt_cluster_expansion_vs_a2_L$(L).pdf"
plot_a2, error_plot_a2 = plot_free_energy_vs_a2_with_single_site_clusters_efficient(
            L, a1, a2_range, max_weights; save_path=save_path_a2)



a1 = 0.5; a2 = 0.5
save_path_error = "../visualization/aklt_error_vs_weight_L$(L)_a1$(a1)_a2$(a2).pdf"
error_plot = plot_error_vs_weight(L, a1, a2, max_weights; 
                                        save_path=save_path_error,
                                        fontfamily="Times New Roman")


a1 = 0.5; a2 = 1.5
save_path_error = "../visualization/aklt_error_vs_weight_L$(L)_a1$(a1)_a2$(a2).pdf"
error_plot = plot_error_vs_weight(L, a1, a2, max_weights; 
                                        save_path=save_path_error,
                                        fontfamily="Times New Roman")


a1 = 0.5; a2 = 5.0
save_path_error = "../visualization/aklt_error_vs_weight_L$(L)_a1$(a1)_a2$(a2).pdf"
error_plot = plot_error_vs_weight(L, a1, a2, max_weights; 
                                        save_path=save_path_error,
                                        fontfamily="Times New Roman")

