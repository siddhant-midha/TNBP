#!/usr/bin/env julia

"""
cluster_expansion_free_energy.jl

Implementation of cluster expansion for 2D Ising model free energy density using saved clusters.
Follows exactly the workflow from functions/cluster_expansion_2d_ising.jl but uses pre-computed clusters.
"""

include("dependencies.jl")
include("functions/ClusterEnumeration.jl")
include("functions/Ising2D.jl")
include("functions/BP.jl")
using Serialization
using Plots
using ITensors, Graphs

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

# Ursell function for connected clusters (from original implementation)
function ursell_function(cluster::Cluster, all_loops::Vector{Loop})
    """
    Compute the Ursell function φ(W) for a connected cluster W.
    For connected clusters, this is (-1)^(|W|-1) * (|W|-1)!
    where |W| is the total number of loops in the cluster (with multiplicities).
    """
    total_loops = cluster.total_loops
    if total_loops == 1
        return 1.0
    else
        # φ(W) = (-1)^(|W|-1) * (|W|-1)!
        sign = (total_loops - 1) % 2 == 0 ? 1.0 : -1.0
        factorial_part = factorial(big(total_loops - 1))
        return sign * Float64(factorial_part)
    end
end

function cluster_expansion_2d_ising_with_saved_clusters(cluster_file::String, L::Int, β::Float64, h::Float64=0.0; 
                                                       max_weights=[4, 6, 8, 10])
    """
    Complete cluster expansion workflow for 2D Ising model using saved clusters.
    Follows exactly the pattern from functions/cluster_expansion_2d_ising.jl
    """
    
    println("="^80)
    println("🔥 Cluster Expansion for 2D Ising Model (Using Saved Clusters)")
    println("="^80)
    println("Parameters: L=$L, β=$β, h=$h")
    println("Cluster weights: $max_weights")
    println()
    
    # Step 1: Create Ising tensor network
    println("Step 1: Creating Ising tensor network...")
    T = Ising2D.get_ising_tn(L, β; h=h)
    N = L^2
    println("✅ Created $N tensors for $(L)×$(L) lattice")
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    println("✅ Found $(length(edges)) edges")
    
    # Step 2: Compute BP fixed point
    println("\nStep 2: Computing BP fixed point...")
    messages = BP.get_messages(T, edges, links)
    messages = BP.message_passing(T, messages, edges, adj_mat; max_iters=1000)
    println("✅ BP converged")
    
    # Get BP partition function (before normalization)
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Z_bp_full)
    println("✅ BP log partition function: $log_Z_bp_full")
    
    # Step 3: Normalize tensors
    println("\nStep 3: Normalizing tensors...")
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)
    
    # Normalization factor
    normalization_factor = sum(log.(real.(Z_l)))
    println("✅ Normalization factor: $normalization_factor")
    
    # Verify normalization (BP on normalized tensors should give 1)
    Z_bp_normalized = BP.mean_free_partition_fn(1:N, T_normalized, messages, adj_mat)
    println("✅ Normalized BP partition function: $Z_bp_normalized (should ≈ 1)")
    
    # Step 4: Load saved clusters and compute corrections
    println("\nStep 4: Loading saved clusters...")
    cluster_data, cluster_summary = load_cluster_data(cluster_file)
    
    results = Dict{String, Any}()
    results["L"] = L
    results["beta"] = β
    results["h"] = h
    results["bp_log_Z"] = log_Z_bp_full  # Keep as raw log partition function
    results["normalization_factor"] = normalization_factor
    results["cluster_corrections"] = Dict()
    
    for max_weight in max_weights
        println("\n" * "="^60)
        println("🔍 Computing cluster expansion up to weight $max_weight")
        println("="^60)
        
        correction = compute_cluster_correction_with_saved_data(
            T_normalized, messages, edges, links, adj_mat, 
            cluster_data, L, max_weight)
        
        # Total log partition function with cluster correction
        # According to polymer theory: log(Z) = log(Z_BP) + Σ φ(W) Z_W
        # For free energy density: f = -log(Z)/(2*N), so correction is added directly
        log_Z_corrected = log_Z_bp_full + correction * N  # correction is per-site, multiply by N for total log(Z)
        
        results["cluster_corrections"][max_weight] = Dict(
            "correction" => correction,  # This is the free energy density correction
            "total_log_Z" => log_Z_corrected,
            "improvement" => correction
        )
        
        println("📊 Results for weight $max_weight:")
        println("  Cluster correction: $correction")
        println("  Total log partition function: $log_Z_corrected")
    end
    
    # Step 5: Compare with exact Onsager solution
    if h == 0.0  # Onsager solution only valid for h=0
        println("\n" * "="^60)
        println("📐 Comparing with exact Onsager solution")
        println("="^60)
        
        exact_free_energy = Ising2D.free_energy(β)
        exact_log_Z = -exact_free_energy * (2 * N)  # Convert free energy back to log Z
        results["exact_onsager"] = exact_free_energy
        results["exact_log_Z"] = exact_log_Z
        
        println("🎯 Exact Onsager free energy: $exact_free_energy")
        println("🎯 Exact log partition function: $exact_log_Z")
        println("\n📈 Comparison with exact result:")
        println("  BP approximation error: $(abs(results["bp_log_Z"] - exact_log_Z))")
        
        for max_weight in max_weights
            corrected_log_Z = results["cluster_corrections"][max_weight]["total_log_Z"]
            error = abs(corrected_log_Z - exact_log_Z)
            improvement = abs(results["bp_log_Z"] - exact_log_Z) - error
            
            println("  Weight $max_weight error: $error (improvement: $improvement)")
        end
    else
        println("⚠️  Exact comparison only available for h=0 (Onsager solution)")
        results["exact_onsager"] = nothing
        results["exact_log_Z"] = nothing
    end
    
    return results
end

function compute_cluster_correction_with_saved_data(T_normalized, messages, edges, links, adj_mat, 
                                                   cluster_data, L::Int, max_weight::Int)
    """
    Compute the cluster expansion correction using saved cluster data.
    
    CRITICAL: The saved clusters are "representatives" from one site with translational 
    symmetry removed. For periodic boundary conditions, we need to multiply by the 
    number of equivalent translations (N = L²) to get the total contribution.
    """
    
    # Get clusters from saved data
    if hasfield(typeof(cluster_data), :clusters_by_site)
        # All-sites enumeration data - use site 1 clusters (PBC makes all equivalent)
        representative_clusters = cluster_data.clusters_by_site[1]
    else
        # Single-site enumeration data
        representative_clusters = cluster_data.clusters
    end
    
    all_loops = cluster_data.all_loops
    N = L^2  # Total number of sites
    
    println("📊 Using $(length(representative_clusters)) representative clusters from saved data")
    println("   Available loops: $(length(all_loops))")
    println("   Will multiply by N=$N for translational equivalence in PBC")
    
    # Filter clusters by weight
    relevant_clusters = [c for c in representative_clusters if c.weight <= max_weight]
    println("   Using $(length(relevant_clusters)) representative clusters (weight ≤ $max_weight)")
    
    if isempty(relevant_clusters)
        println("⚠️  No relevant clusters found, returning zero correction")
        return 0.0
    end
    
    # Compute total correction from representatives
    representative_correction = 0.0
    
    println("💫 Computing cluster contributions...")
    # Process all relevant clusters (no arbitrary limit like in demo)
    
    for (i, cluster) in enumerate(relevant_clusters)
        if i % 100 == 0 || i <= 10
            println("  Processing representative cluster $i/$(length(relevant_clusters))...")
        end
        
        # Compute Ursell function φ(W)
        phi_W = ursell_function(cluster, all_loops)
        
        if abs(phi_W) < 1e-15
            continue  # Skip negligible contributions
        end
        
        # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
        Z_W = 1.0 + 0im  # Complex number for cluster contribution
        
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
                continue
            end
        end
        
        # Add weighted contribution: φ(W) * Z_W
        contribution = phi_W * Z_W
        representative_correction += contribution
        
        if abs(contribution) > 1e-10 && i <= 20
            println("    Representative cluster $i: φ=$phi_W, Z_W=$Z_W, contribution=$contribution")
        end
    end
    
    # CRITICAL: For PBC, multiply by N to account for all translational copies
    # But since we're computing FREE ENERGY DENSITY, this factor of N cancels out
    # when we divide by N in the final free energy density calculation.
    # So the total correction to the LOG partition function is N * representative_correction,
    # but the correction to the FREE ENERGY DENSITY is just representative_correction.
    total_log_Z_correction = N * representative_correction
    free_energy_density_correction = representative_correction  # This is what gets returned
    
    println("✅ Representative cluster correction: $representative_correction")
    println("✅ Total log(Z) correction: $total_log_Z_correction")
    println("✅ Free energy density correction: $free_energy_density_correction")
    return real(free_energy_density_correction)  # Return free energy density correction
end

# Plot free energy vs beta
function plot_free_energy_vs_beta_with_saved_clusters(cluster_file::String, β_range::Vector{Float64}, 
                                                     max_weights::Vector{Int}=[4,6,8,10]; save_path::String="")
    """
    Plot free energy density as a function of β using saved clusters.
    """
    
    println("📊 Creating free energy vs β plot using saved clusters")
    println("   β range: $(β_range[1]) to $(β_range[end]) ($(length(β_range)) points)")
    println("   Weight truncations: $max_weights")
    
    # Load cluster data to get lattice size
    cluster_data, cluster_summary = load_cluster_data(cluster_file)
    L = cluster_summary["lattice_size"]
    
    # Initialize plot
    p = plot(title="2D Ising Free Energy Density vs β (Cluster Expansion)", 
             xlabel="β", ylabel="Free Energy Density f",
             legend=:topleft, size=(800, 600),
             dpi=300)
    
    # Compute exact Onsager solution
    println("\n🎯 Computing exact Onsager solution...")
    exact_f = [Ising2D.free_energy(β) for β in β_range]
    plot!(p, β_range, exact_f, linewidth=3, color=:black, 
          label="Exact (Onsager)", linestyle=:solid)
    
    # Compute BP baseline (more efficiently - just compute BP without clusters)
    println("\n📊 Computing BP baseline...")
    bp_f = Float64[]
    for (j, β) in enumerate(β_range)
        if j % 5 == 1
            println("   BP β = $β ($(j)/$(length(β_range)))")
        end
        
        try
            # Quick BP computation without clusters
            T = Ising2D.get_ising_tn(L, β; h=0.0)
            N = L^2
            adj_mat, edges, links = BP.get_adj_mat(T)
            messages = BP.get_messages(T, edges, links)
            messages = BP.message_passing(T, messages, edges, adj_mat; max_iters=1000)
            Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
            
            # BP free energy density: f = -log(Z_BP) / (2*N)
            bp_log_Z = log(Z_bp_full)
            f = -real(bp_log_Z) / (2 * L^2)  # Divide by N for BP
            push!(bp_f, f)
        catch e
            println("     ⚠️ Error at β=$β: $e")
            push!(bp_f, NaN)
        end
    end
    
    # Plot BP baseline
    plot!(p, β_range, bp_f, linewidth=2, color=:gray, 
          label="BP (no corrections)", linestyle=:dot)
    
    # Colors for different weights
    colors = [:red, :blue, :green, :orange, :purple, :brown]
    
    # Compute cluster expansion for each weight
    for (i, max_weight) in enumerate(max_weights)
        println("\n📈 Computing cluster expansion for weight $max_weight...")
        
        cluster_f = Float64[]
        
        for (j, β) in enumerate(β_range)
            if j % 5 == 1
                println("   β = $β ($(j)/$(length(β_range)))")
            end
            
            try
                results = cluster_expansion_2d_ising_with_saved_clusters(cluster_file, L, β, 0.0; max_weights=[max_weight])
                # For cluster expansion: f = -log(Z_BP) / (2*N) + correction
                # where correction is already the free energy density correction
                bp_log_Z = results["bp_log_Z"]
                correction = results["cluster_corrections"][max_weight]["correction"]
                # Both BP and correction are now free energy density terms
                f = -real(bp_log_Z) / (2 * L^2) + real(correction)
                push!(cluster_f, f)
            catch e
                println("     ⚠️ Error at β=$β: $e")
                push!(cluster_f, NaN)
            end
        end
        
        # Plot this weight truncation
        color = colors[min(i, length(colors))]
        plot!(p, β_range, cluster_f, linewidth=2, color=color,
              label="Weight ≤ $max_weight", linestyle=:dash)
    end
    
    # Save plot
    if !isempty(save_path)
        println("\n💾 Saving plot to: $save_path")
        savefig(p, save_path)
    end
    
    return p
end

# Main execution function
function main()
    println("🚀 Cluster Expansion Free Energy Calculator (Using Saved Clusters)")
    println("="^80)
    
    # Parameters
    cluster_file = "saved_clusters/periodic_clusters_L11_w10_2025-08-14T20-25-54-769.jld2"
    β_range = collect(0.1:0.2:1.6)  # β from 0.1 to 1.6 with fewer points for faster computation
    max_weights = [4, 6]  # Start with smaller weights to ensure completion
    
    println("📋 Configuration:")
    println("   Cluster file: $(basename(cluster_file))")
    println("   β range: $(β_range[1]) to $(β_range[end]) ($(length(β_range)) points)")
    println("   Weight truncations: $max_weights")
    
    # Create plot
    save_path = "visualization/cluster_expansion_free_energy_L11_corrected.png"
    plot_obj = plot_free_energy_vs_beta_with_saved_clusters(cluster_file, β_range, max_weights; save_path=save_path)
    
    # Display final results
    println("\n✅ Analysis complete!")
    println("📊 Plot saved to: $save_path")
    
    return plot_obj
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end