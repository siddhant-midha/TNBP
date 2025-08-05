"""
final_comprehensive_plot.jl

Generate the comprehensive free energy vs β plot with all requested cluster weights.
"""

include("efficient_beta_plot.jl")

function generate_final_plot()
    """Generate the final comprehensive plot as requested."""
    
    println("🎯 Generating Final Comprehensive Free Energy vs β Plot")
    println("="^80)
    
    # Parameters for good quality plot
    L = 5  # System size (compromise between accuracy and computational cost)
    n_points = 12  # Good resolution across β range
    max_weights = [4, 6, 8, 10, 12]  # All requested weights
    
    # Run the comprehensive analysis
    results = efficient_beta_plot(L, n_points; max_weights=max_weights)
    
    # Create additional focused plot around critical region
    beta_critical = range(0.35, 0.55, length=8)
    
    println("\n🔍 Creating focused critical region plot...")
    
    # Extract data for focused plot
    critical_bp = Float64[]
    critical_exact = Float64[]
    critical_clusters = Dict{Int, Vector{Float64}}()
    
    for weight in max_weights
        critical_clusters[weight] = Float64[]
    end
    
    for β in beta_critical
        println("  Critical region β = $β")
        try
            # Quick computation for critical region
            T = Ising2D.get_ising_tn(L, β; h=0.0)
            adj_mat, edges, links = BP.get_adj_mat(T)
            messages = BP.get_messages(T, edges, links)
            messages = BP.message_passing(T, messages, edges, links, adj_mat; max_iters=1000)
            
            Z_bp = BP.mean_free_partition_fn(1:L^2, T, messages, edges, links, adj_mat)
            bp_f = -log(Z_bp) / (2 * L^2)
            push!(critical_bp, real(bp_f))
            
            exact_f = Ising2D.free_energy(β)
            push!(critical_exact, exact_f)
            
            Z_l = BP.get_fixed_point_list(T, messages, edges, links, adj_mat)
            T_norm = BP.normalize_tensor(T, Z_l)
            
            for weight in max_weights
                try
                    if weight <= 100
                        correction = compute_limited_cluster_correction(T_norm, messages, edges, links, adj_mat, L, weight)
                        corrected_f = -(log(Z_bp) + correction) / (2 * L^2)
                        push!(critical_clusters[weight], real(corrected_f))
                    else
                        # Extrapolate for higher weights
                        base_correction = critical_clusters[6][end] - bp_f
                        estimated_f = bp_f + base_correction * (weight/6)^0.3
                        push!(critical_clusters[weight], estimated_f)
                    end
                catch e
                    push!(critical_clusters[weight], NaN)
                end
            end
        catch e
            println("    Error: $e")
            push!(critical_bp, NaN)
            push!(critical_exact, NaN)
            for weight in max_weights
                push!(critical_clusters[weight], NaN)
            end
        end
    end
    
    # Create focused critical region plot
    p_critical = plot(title="Free Energy Near Critical Point (β_c ≈ 0.44)",
                     xlabel="β (Inverse Temperature)", 
                     ylabel="Free Energy Density",
                     legend=:bottomleft,
                     size=(800, 600),
                     dpi=300)
    
    plot!(p_critical, beta_critical, critical_exact, 
          label="Exact (Onsager)", linewidth=4, color=:black)
    
    plot!(p_critical, beta_critical, critical_bp,
          label="BP Approximation", linewidth=3, color=:red, linestyle=:dash)
    
    colors = [:blue, :green, :orange, :purple, :brown]
    for (i, weight) in enumerate(max_weights)
        if i <= length(colors)
            plot!(p_critical, beta_critical, critical_clusters[weight],
                  label="Weight $weight", linewidth=2, color=colors[i])
        end
    end
    
    # Add critical point marker
    β_c = log(1 + sqrt(2))/2
    plot!(p_critical, [β_c, β_c], [-1.2, -0.6], 
          label="β_c (exact)", color=:gray, linestyle=:dot, linewidth=2)
    
    # Save critical region plot
    savefig(p_critical, "critical_region_free_energy.png")
    println("✅ Critical region plot saved as: critical_region_free_energy.png")
    
    println("\n📈 FINAL SUMMARY")
    println("="^50)
    println("System analyzed: $(L)x$(L) 2D Ising model")
    println("β range: 0.1 to 2.0")
    println("Cluster weights: $max_weights")
    println("\nNear critical point (β ≈ $β_c):")
    
    # Find closest point to critical
    critical_idx = argmin(abs.(beta_critical .- β_c))
    β_actual = beta_critical[critical_idx]
    
    println("  β = $β_actual")
    println("  Exact free energy: $(critical_exact[critical_idx])")
    println("  BP approximation: $(critical_bp[critical_idx])")
    bp_error = abs(critical_bp[critical_idx] - critical_exact[critical_idx])
    println("  BP error: $bp_error")
    
    for weight in max_weights
        if !isnan(critical_clusters[weight][critical_idx])
            cluster_f = critical_clusters[weight][critical_idx]
            cluster_error = abs(cluster_f - critical_exact[critical_idx])
            improvement = bp_error - cluster_error
            println("  Weight $weight: $cluster_f (error: $cluster_error, improvement: $improvement)")
        end
    end
    
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    final_results = generate_final_plot()
end