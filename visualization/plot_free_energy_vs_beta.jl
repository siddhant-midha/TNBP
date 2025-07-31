"""
plot_free_energy_vs_beta.jl

Generate plots of free energy density as a function of β (inverse temperature)
comparing cluster expansion corrections at different weights with exact Onsager solution.
"""

include("../functions/cluster_expansion_2d_ising.jl")
using Plots

function plot_free_energy_vs_beta(L::Int=6, beta_range=(0.1, 2.0), n_points::Int=20;
                                  max_weights=[4, 6, 8, 10, 12])
    """
    Plot free energy density vs β for different cluster expansion weights.
    
    Args:
        L: Lattice size (default: 6x6 for computational efficiency)
        beta_range: Range of β values (default: 0.1 to 2.0)
        n_points: Number of β points to compute (default: 20)
        max_weights: Cluster weights to compare (default: [4,6,8,10,12])
    """
    
    println("="^80)
    println("🔥 Plotting Free Energy vs β for 2D Ising Model")
    println("="^80)
    println("System: $(L)x$(L) lattice")
    println("β range: $(beta_range[1]) to $(beta_range[2])")
    println("Points: $n_points")
    println("Cluster weights: $max_weights")
    println()
    
    # Create β values
    beta_values = range(beta_range[1], beta_range[2], length=n_points)
    
    # Initialize storage arrays
    bp_free_energies = Float64[]
    exact_free_energies = Float64[]
    cluster_corrections = Dict{Int, Vector{Float64}}()
    
    for weight in max_weights
        cluster_corrections[weight] = Float64[]
    end
    
    # Compute free energies for each β
    println("Computing free energies...")
    for (i, β) in enumerate(beta_values)
        println("Progress: $i/$n_points (β = $β)")
        
        try
            # Run cluster expansion for this β
            results = cluster_expansion_2d_ising(L, β, 0.0; max_weights=max_weights)
            
            # Store results
            push!(bp_free_energies, results["bp_log_Z"])
            push!(exact_free_energies, results["exact_onsager"])
            
            for weight in max_weights
                if haskey(results["cluster_corrections"], weight)
                    push!(cluster_corrections[weight], 
                          results["cluster_corrections"][weight]["total_free_energy"])
                else
                    push!(cluster_corrections[weight], NaN)
                end
            end
            
        catch e
            println("⚠️  Error at β=$β: $e")
            # Fill with NaN for failed computations
            push!(bp_free_energies, NaN)
            push!(exact_free_energies, NaN)
            for weight in max_weights
                push!(cluster_corrections[weight], NaN)
            end
        end
    end
    
    # Create plots
    println("\n📊 Creating plots...")
    
    # Main comparison plot
    p1 = plot(title="Free Energy Density vs Inverse Temperature",
              xlabel="β (Inverse Temperature)", 
              ylabel="Free Energy Density",
              legend=:topright,
              size=(800, 600),
              dpi=300)
    
    # Plot exact solution
    plot!(p1, beta_values, exact_free_energies, 
          label="Exact (Onsager)", linewidth=3, color=:black, linestyle=:solid)
    
    # Plot BP approximation
    plot!(p1, beta_values, bp_free_energies,
          label="BP Approximation", linewidth=2, color=:red, linestyle=:dash)
    
    # Plot cluster corrections
    colors = [:blue, :green, :orange, :purple, :brown]
    for (i, weight) in enumerate(max_weights)
        if i <= length(colors)
            plot!(p1, beta_values, cluster_corrections[weight],
                  label="Cluster Weight $weight", linewidth=2, 
                  color=colors[i], linestyle=:solid)
        end
    end
    
    # Error plot
    p2 = plot(title="Absolute Error vs Exact Solution",
              xlabel="β (Inverse Temperature)", 
              ylabel="|Free Energy - Exact|",
              legend=:topright,
              size=(800, 600),
              yscale=:log10,
              dpi=300)
    
    # Compute and plot errors
    bp_errors = abs.(bp_free_energies .- exact_free_energies)
    plot!(p2, beta_values, bp_errors,
          label="BP Error", linewidth=2, color=:red, linestyle=:dash)
    
    for (i, weight) in enumerate(max_weights)
        if i <= length(colors)
            cluster_errors = abs.(cluster_corrections[weight] .- exact_free_energies)
            plot!(p2, beta_values, cluster_errors,
                  label="Weight $weight Error", linewidth=2, 
                  color=colors[i], linestyle=:solid)
        end
    end
    
    # Improvement plot
    p3 = plot(title="Improvement Over BP",
              xlabel="β (Inverse Temperature)", 
              ylabel="Error Reduction",
              legend=:topright,
              size=(800, 600),
              dpi=300)
    
    for (i, weight) in enumerate(max_weights)
        if i <= length(colors)
            cluster_errors = abs.(cluster_corrections[weight] .- exact_free_energies)
            improvements = bp_errors .- cluster_errors
            plot!(p3, beta_values, improvements,
                  label="Weight $weight Improvement", linewidth=2, 
                  color=colors[i], linestyle=:solid)
        end
    end
    
    # Combine plots
    combined_plot = plot(p1, p2, p3, layout=(3,1), size=(800, 1200))
    
    # Save plots
    output_file = "free_energy_vs_beta_L$(L).png"
    savefig(combined_plot, output_file)
    println("✅ Plot saved as: $output_file")
    
    # Also save individual plots
    savefig(p1, "free_energy_comparison_L$(L).png")
    savefig(p2, "error_comparison_L$(L).png") 
    savefig(p3, "improvement_comparison_L$(L).png")
    
    # Print summary statistics
    println("\n📈 Summary Statistics:")
    println("β range: $(minimum(beta_values)) to $(maximum(beta_values))")
    
    # Find critical region (around β_c ≈ 0.44 for 2D Ising)
    critical_idx = argmin(abs.(beta_values .- 0.44))
    β_c = beta_values[critical_idx]
    
    println("Near critical point (β ≈ $β_c):")
    println("  Exact free energy: $(exact_free_energies[critical_idx])")
    println("  BP error: $(bp_errors[critical_idx])")
    
    for weight in max_weights
        if !isnan(cluster_corrections[weight][critical_idx])
            error = abs(cluster_corrections[weight][critical_idx] - exact_free_energies[critical_idx])
            improvement = bp_errors[critical_idx] - error
            println("  Weight $weight error: $error (improvement: $improvement)")
        end
    end
    
    return Dict(
        "beta_values" => beta_values,
        "bp_free_energies" => bp_free_energies,
        "exact_free_energies" => exact_free_energies,
        "cluster_corrections" => cluster_corrections,
        "plot" => combined_plot
    )
end

function quick_beta_scan(L::Int=4, n_points::Int=10)
    """Quick beta scan with smaller system for testing."""
    
    println("🚀 Quick β scan with $(L)x$(L) system")
    return plot_free_energy_vs_beta(L, (0.2, 1.0), n_points; max_weights=[4, 6])
end

function detailed_beta_scan(L::Int=6, n_points::Int=15)
    """Detailed beta scan with moderate system size."""
    
    println("🔬 Detailed β scan with $(L)x$(L) system")
    return plot_free_energy_vs_beta(L, (0.1, 2.0), n_points; max_weights=[4, 6, 8, 10])
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("🎯 Generating free energy vs β plots...")
    
    # Run quick test first
    println("Running quick test...")
    quick_results = quick_beta_scan(4, 8)
    
    println("\nRunning detailed analysis...")
    detailed_results = detailed_beta_scan(6, 12)
    
    println("\n✅ All plots generated successfully!")
    println("Check the generated .png files for results.")
end