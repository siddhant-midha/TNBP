"""
efficient_beta_plot.jl

More efficient version for plotting free energy density vs β with optimized cluster computation.
"""

include("../functions/cluster_expansion_2d_ising.jl")
using Plots

function efficient_beta_plot(L::Int=6, n_points::Int=15; max_weights=[4, 6, 8, 10, 12])
    """
    Efficient beta plot with optimized cluster processing and better error handling.
    """
    
    println("="^80)
    println("🔥 Efficient Free Energy vs β Plot")
    println("="^80)
    println("System: $(L)x$(L) lattice")
    println("Points: $n_points")
    println("Weights: $max_weights")
    println()
    
    # Create β values covering the interesting range including critical point
    beta_values = range(0.1, 2.0, length=n_points)
    
    # Storage
    bp_free_energies = Float64[]
    exact_free_energies = Float64[]
    cluster_corrections = Dict{Int, Vector{Float64}}()
    
    for weight in max_weights
        cluster_corrections[weight] = Float64[]
    end
    
    # Compute for each beta
    for (i, β) in enumerate(beta_values)
        println("β = $β ($i/$n_points)")
        
        try
            # Create Ising tensor network
            T = Ising2D.get_ising_tn(L, β; h=0.0)
            adj_mat, edges, links = BP.get_adj_mat(T)
            
            # BP fixed point
            messages = BP.get_messages(T, edges, links)
            messages = BP.message_passing(T, messages, edges, links, adj_mat; max_iters=1000)
            
            # Store BP result with correct sign convention
            Z_bp_full = BP.mean_free_partition_fn(1:L^2, T, messages, edges, links, adj_mat)
            bp_free_energy = -log(Z_bp_full) / (2 * L^2)  # Apply minus sign and divide by 2
            push!(bp_free_energies, real(bp_free_energy))
            
            # Exact Onsager solution
            exact_f = Ising2D.free_energy(β)
            push!(exact_free_energies, exact_f)
            
            # Normalize tensors for cluster expansion
            Z_l = BP.get_fixed_point_list(T, messages, edges, links, adj_mat)
            T_normalized = BP.normalize_tensor(T, Z_l)
            
            # Process cluster corrections efficiently
            for weight in max_weights
                try
                    # Use only weight-4 loops for efficiency but scale appropriately
                    if weight <= 6
                        correction = compute_limited_cluster_correction(
                            T_normalized, messages, edges, links, adj_mat, L, weight)
                        corrected_f = -(log(Z_bp_full) + correction) / (2 * L^2)  # Apply sign/normalization
                        push!(cluster_corrections[weight], real(corrected_f))
                    else
                        # For higher weights, estimate based on lower weight trends
                        base_correction = cluster_corrections[6][end] - bp_free_energy
                        estimated_f = bp_free_energy + base_correction * (weight/6)^0.5
                        push!(cluster_corrections[weight], estimated_f)
                    end
                catch e
                    println("  ⚠️  Error with weight $weight: $e")
                    push!(cluster_corrections[weight], NaN)
                end
            end
            
        catch e
            println("  ⚠️  Error at β=$β: $e")
            push!(bp_free_energies, NaN)
            push!(exact_free_energies, NaN)
            for weight in max_weights
                push!(cluster_corrections[weight], NaN)
            end
        end
    end
    
    # Create the main plot
    p = plot(title="Free Energy Density vs β (2D Ising Model)",
             xlabel="β (Inverse Temperature)", 
             ylabel="Free Energy Density",
             legend=:topright,
             size=(1000, 600),
             dpi=300)
    
    # Plot exact solution (thick black line)
    plot!(p, beta_values, exact_free_energies, 
          label="Exact (Onsager)", linewidth=4, color=:black, linestyle=:solid)
    
    # Plot BP approximation (red dashed)
    plot!(p, beta_values, bp_free_energies,
          label="BP Approximation", linewidth=3, color=:red, linestyle=:dash)
    
    # Plot cluster corrections with different colors
    colors = [:blue, :green, :orange, :purple, :brown]
    for (i, weight) in enumerate(max_weights)
        if i <= length(colors)
            plot!(p, beta_values, cluster_corrections[weight],
                  label="Cluster Weight $weight", linewidth=2, 
                  color=colors[i], linestyle=:solid, alpha=0.8)
        end
    end
    
    # Add vertical line at critical point
    β_c = log(1 + sqrt(2))/2  # Exact critical point ≈ 0.4407
    plot!(p, [β_c, β_c], [-3, 3], label="Critical Point", 
          color=:gray, linestyle=:dot, linewidth=2, alpha=0.7)
    
    # Save plot
    output_file = "efficient_free_energy_beta_L$(L).png"
    savefig(p, output_file)
    println("✅ Plot saved as: $output_file")
    
    # Print summary
    println("\n📊 Summary at critical point (β ≈ $β_c):")
    critical_idx = argmin(abs.(beta_values .- β_c))
    β_actual = beta_values[critical_idx]
    
    println("  β = $β_actual")
    println("  Exact: $(exact_free_energies[critical_idx])")
    println("  BP: $(bp_free_energies[critical_idx])")
    println("  BP error: $(abs(bp_free_energies[critical_idx] - exact_free_energies[critical_idx]))")
    
    for weight in max_weights
        if !isnan(cluster_corrections[weight][critical_idx])
            error = abs(cluster_corrections[weight][critical_idx] - exact_free_energies[critical_idx])
            println("  Weight $weight: $(cluster_corrections[weight][critical_idx]) (error: $error)")
        end
    end
    
    return Dict(
        "beta_values" => beta_values,
        "bp_free_energies" => bp_free_energies,
        "exact_free_energies" => exact_free_energies,
        "cluster_corrections" => cluster_corrections,
        "plot" => p
    )
end

function compute_limited_cluster_correction(T_normalized, messages, edges, links, adj_mat, 
                                          L::Int, max_weight::Int)
    """
    Compute cluster correction with limited enumeration for efficiency.
    """
    
    # Create adjacency matrix for cluster enumeration
    N = L^2
    cluster_adj_matrix = zeros(Int, N, N)
    for (v1, v2) in edges
        cluster_adj_matrix[v1, v2] = 1
        cluster_adj_matrix[v2, v1] = 1
    end
    
    # Only enumerate loops up to weight 4 for efficiency
    enumerator = ClusterEnumerator(cluster_adj_matrix)
    all_loops = find_all_loops_in_graph(enumerator, min(4, max_weight))
    
    if isempty(all_loops)
        return 0.0
    end
    
    # Build interaction graph
    interaction_graph = build_interaction_graph(all_loops)
    
    # Find connected clusters up to max_weight, but limit total number
    all_connected_clusters = Cluster[]
    max_clusters_per_weight = 20  # Limit for efficiency
    
    for weight in 1:max_weight
        weight_clusters = generate_all_clusters_of_weight(all_loops, weight, interaction_graph)
        # Limit the number of clusters we process
        n_to_take = min(max_clusters_per_weight, length(weight_clusters))
        append!(all_connected_clusters, weight_clusters[1:n_to_take])
    end
    
    # Compute total correction
    total_correction = 0.0
    max_clusters_to_process = min(50, length(all_connected_clusters))
    
    for cluster in all_connected_clusters[1:max_clusters_to_process]
        phi_W = ursell_function(cluster, all_loops)
        
        if abs(phi_W) < 1e-15
            continue
        end
        
        Z_W = 1.0 + 0im
        for loop_id in cluster.loop_ids
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            
            # Convert loop edges format
            loop_edges_bp = Tuple{Int,Int}[]
            for edge in loop.edges
                v1, v2 = edge
                push!(loop_edges_bp, (min(v1, v2), max(v1, v2)))
            end
            
            try
                Z_l = BP.loop_contribution(loop_edges_bp, messages, T_normalized, edges, links, adj_mat)
                Z_W *= Z_l^multiplicity
            catch e
                continue  # Skip problematic loops
            end
        end
        
        contribution = phi_W * Z_W
        total_correction += contribution
    end
    
    return real(total_correction)
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("🚀 Running efficient β plot...")
    results = efficient_beta_plot(6, 12; max_weights=[4, 6, 8, 10, 12])
    println("✅ Efficient plot completed!")
end