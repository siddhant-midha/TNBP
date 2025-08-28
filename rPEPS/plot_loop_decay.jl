using JLD2, FileIO, Plots, Statistics
using Printf

"""
Plotting script for loop decay analysis.
Creates plots showing how loop contributions decay with loop weight for different η and N values.
"""

function load_loop_decay_results(results_dir::String="loop_decay_results")
    """
    Load and organize loop decay results from JLD2 files.
    """
    if !isdir(results_dir)
        println("Warning: Directory $results_dir not found!")
        return Dict()
    end
    
    files = filter(f -> endswith(f, ".jld2"), readdir(results_dir))
    results = Dict()
    
    println("Found $(length(files)) loop decay result files")
    
    for file in files
        try
            filepath = joinpath(results_dir, file)
            data = load(filepath, "results")
            N, η = data["N"], data["eta"]
            
            if !haskey(results, N)
                results[N] = Dict()
            end
            results[N][η] = data
            
        catch e
            println("Warning: Failed to load $file: $e")
        end
    end
    
    return results
end

function plot_loop_decay_by_N(results::Dict)
    """
    Create separate plots for each N showing loop decay vs weight for different η values.
    """
    for N in sort(collect(keys(results)))
        println("Creating plot for N=$N...")
        
        # Get all η values for this N
        η_values = sort(collect(keys(results[N])))
        
        if isempty(η_values)
            println("No η values for N=$N, skipping...")
            continue
        end
        
        # Create plot
        p = plot(title="Loop Contribution Decay vs Weight (N=$N)", 
                 xlabel="Loop Weight", 
                 ylabel="Average |Loop Contribution|",
                 yscale=:log10,
                 legend=:topright, 
                 size=(900, 600), 
                 dpi=300)
        
        # Color palette for different η values
        colors = palette(:viridis, length(η_values))
        
        for (i, η) in enumerate(η_values)
            data = results[N][η]
            avg_contributions = data["avg_contributions_by_weight"]
            
            if !isempty(avg_contributions)
                weights = sort(collect(keys(avg_contributions)))
                # Filter to only even weights
                even_weights = filter(w -> w % 2 == 0, weights)
                contributions = [avg_contributions[w] for w in even_weights]
                
                if !isempty(even_weights)  # Only plot if we have even weights
                    # Add small offset to avoid log(0)
                    contributions = contributions .+ 1e-50
                    
                    plot!(p, even_weights, contributions,
                          linewidth=2,
                          marker=:circle,
                          markersize=4,
                          label="η = $η",
                          color=colors[i])
                end
            end
        end
        
        # Save plot
        figs_dir = ensure_figs_dir()
        filename = joinpath(figs_dir, "loop_decay_N$(N).png")
        savefig(p, filename)
        display(p)
        println("Saved: $filename")
    end
end

function plot_loop_decay_by_eta(results::Dict, target_etas::Vector{Float64}=[0.2, 0.5, 0.8, 1.2, 1.8])
    """
    Create plots for specific η values showing how decay varies with N.
    """
    for target_η in target_etas
        println("Creating plot for η≈$target_η...")
        
        p = plot(title="Loop Contribution Decay vs Weight (η≈$target_η)", 
                 xlabel="Loop Weight", 
                 ylabel="Average |Loop Contribution|",
                 yscale=:log10,
                 legend=:topright, 
                 size=(900, 600), 
                 dpi=300)
        
        colors = [:red, :blue, :green, :purple, :orange]
        markers = [:circle, :square, :diamond, :utriangle, :star5]
        
        plot_created = false
        
        for (i, N) in enumerate(sort(collect(keys(results))))
            # Find closest η to target
            available_etas = collect(keys(results[N]))
            if isempty(available_etas)
                continue
            end
            
            closest_η = available_etas[argmin(abs.(available_etas .- target_η))]
            
            # Only plot if close enough to target
            if abs(closest_η - target_η) <= 0.1
                data = results[N][closest_η]
                avg_contributions = data["avg_contributions_by_weight"]
                
                if !isempty(avg_contributions)
                    weights = sort(collect(keys(avg_contributions)))
                    # Filter to only even weights
                    even_weights = filter(w -> w % 2 == 0, weights)
                    contributions = [avg_contributions[w] for w in even_weights]
                    
                    if !isempty(even_weights)  # Only plot if we have even weights
                        # Add small offset to avoid log(0)
                        contributions = contributions .+ 1e-50
                        
                        color_idx = min(i, length(colors))
                        marker_idx = min(i, length(markers))
                        
                        plot!(p, even_weights, contributions,
                              linewidth=2,
                              marker=markers[marker_idx],
                              markersize=4,
                              label="N = $N (η=$closest_η)",
                              color=colors[color_idx])
                        
                        plot_created = true
                    end
                end
            end
        end
        
        if plot_created
            # Save plot
            figs_dir = ensure_figs_dir()
            filename = joinpath(figs_dir, "loop_decay_eta$(replace(string(target_η), "." => "p")).png")
            savefig(p, filename)
            display(p)
            println("Saved: $filename")
        else
            println("No data found for η≈$target_η")
        end
    end
end

function plot_decay_rate_analysis(results::Dict)
    """
    Analyze and plot the exponential decay rates as a function of η for different N.
    """
    p = plot(title="Loop Decay Rate vs η", 
             xlabel="η", 
             ylabel="Decay Rate (log scale)",
             legend=:topright, 
             size=(900, 600), 
             dpi=300)
    
    colors = [:red, :blue, :green, :purple, :orange]
    markers = [:circle, :square, :diamond, :utriangle, :star5]
    
    for (i, N) in enumerate(sort(collect(keys(results))))
        η_values = []
        decay_rates = []
        
        for η in sort(collect(keys(results[N])))
            data = results[N][η]
            avg_contributions = data["avg_contributions_by_weight"]
            
            if length(avg_contributions) >= 3  # Need at least 3 points for fitting
                weights = sort(collect(keys(avg_contributions)))
                # Filter to only even weights
                even_weights = filter(w -> w % 2 == 0, weights)
                contributions = [avg_contributions[w] for w in even_weights]
                
                # Fit exponential decay: log(contribution) = a - b*weight
                if length(even_weights) >= 3 && all(contributions .> 0)  # Need at least 3 even weights and positive values
                    log_contributions = log.(contributions)
                    
                    # Simple linear fit in log space
                    n = length(even_weights)
                    sum_w = sum(even_weights)
                    sum_logc = sum(log_contributions)
                    sum_w2 = sum(even_weights.^2)
                    sum_w_logc = sum(even_weights .* log_contributions)
                    
                    # Slope (decay rate)
                    decay_rate = -(n * sum_w_logc - sum_w * sum_logc) / (n * sum_w2 - sum_w^2)
                    
                    if decay_rate > 0  # Only plot positive decay rates
                        push!(η_values, η)
                        push!(decay_rates, decay_rate)
                    end
                end
            end
        end
        
        if !isempty(decay_rates)
            color_idx = min(i, length(colors))
            marker_idx = min(i, length(markers))
            
            plot!(p, η_values, decay_rates,
                  linewidth=2,
                  marker=markers[marker_idx],
                  markersize=4,
                  label="N = $N",
                  color=colors[color_idx])
        end
    end
    
    # Save plot
    filename = "loop_decay_rates_vs_eta.png"
    savefig(p, filename)
    display(p)
    println("Saved: $filename")
end

function ensure_figs_dir()
    figs_dir = "figs"
    if !isdir(figs_dir)
        mkpath(figs_dir)
    end
    return figs_dir
end

function main()
    println("="^60)
    println("🔄 Loop Decay Analysis Plotting Script")
    println("="^60)
    
    println("Loading loop decay results...")
    results = load_loop_decay_results()
    
    if isempty(results)
        println("❌ No loop decay results found!")
        println("Make sure loop_decay_results directory exists with .jld2 files.")
        return
    end
    
    println("\n📊 Available data:")
    for N in sort(collect(keys(results)))
        η_vals = sort(collect(keys(results[N])))
        println("  N=$N: η values = $η_vals")
    end
    
    println("\n🎨 Generating plots...")
    
    println("\nGenerating plots by N (different η curves)...")
    plot_loop_decay_by_N(results)
    
    println("\nGenerating plots by η (different N curves)...")
    plot_loop_decay_by_eta(results)
    
    println("\nGenerating decay rate analysis...")
    plot_decay_rate_analysis(results)
    
    println("\n✅ Done! Check the generated PNG files.")
    println("Files generated:")
    println("  - loop_decay_N*.png")
    println("  - loop_decay_eta*.png") 
    println("  - loop_decay_rates_vs_eta.png")
end

# Run the script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
