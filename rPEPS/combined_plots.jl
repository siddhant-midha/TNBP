using JLD2, FileIO, Plots, Statistics
using Printf

"""
Combined plotting script for both cluster and loop correction analysis.
Generates free energy per site error comparison plots showing both methods:
1. BP free energy per site error vs η compared with cluster/loop-corrected errors for different weights
2. Free energy per site error vs weight order for multiple η values showing both cluster and loop corrections
Compatible with old partition function errors and new free energy per site errors.
Cluster corrections: solid lines
Loop corrections: dotted lines
"""

function load_cluster_results(results_dir::String="cluster_correction_results")
    """
    Load and organize cluster correction results from JLD2 files.
    Handles missing files gracefully.
    """
    if !isdir(results_dir)
        println("Warning: Directory $results_dir not found!")
        return Dict()
    end
    
    files = filter(f -> endswith(f, ".jld2"), readdir(results_dir))
    results = Dict()
    
    println("Found $(length(files)) cluster correction result files")
    
    for file in files
        try
            filepath = joinpath(results_dir, file)
            data = load(filepath, "results")
            N, w, η = data["N"], data["w"], data["eta"]
            
            if !haskey(results, N)
                results[N] = Dict()
            end
            if !haskey(results[N], w)
                results[N][w] = Dict()
            end
            results[N][w][η] = data
            
        catch e
            println("Warning: Failed to load cluster file $file: $e")
        end
    end
    
    return results
end

function load_loop_results(results_dir::String="loop_correction_results")
    """
    Load and organize loop correction results from JLD2 files.
    Handles missing files gracefully.
    """
    if !isdir(results_dir)
        println("Warning: Directory $results_dir not found!")
        return Dict()
    end
    
    files = filter(f -> endswith(f, ".jld2"), readdir(results_dir))
    results = Dict()
    
    println("Found $(length(files)) loop correction result files")
    
    for file in files
        try
            filepath = joinpath(results_dir, file)
            data = load(filepath, "results")
            N, w, η = data["N"], data["w"], data["eta"]
            
            if !haskey(results, N)
                results[N] = Dict()
            end
            if !haskey(results[N], w)
                results[N][w] = Dict()
            end
            results[N][w][η] = data
            
        catch e
            println("Warning: Failed to load loop file $file: $e")
        end
    end
    
    return results
end

function get_error_value(data::Dict, error_type::String)
    """
    Extract error value from data with fallback to old keys.
    error_type should be "bp", "cluster", or "loop"
    """
    if error_type == "bp"
        if haskey(data, "bp_fe_error_mean")
            return data["bp_fe_error_mean"]
        elseif haskey(data, "bp_pf_error_mean")
            return data["bp_pf_error_mean"]
        elseif haskey(data, "bp_error_mean")
            return data["bp_error_mean"]
        end
    elseif error_type == "cluster"
        if haskey(data, "cluster_fe_error_mean")
            return data["cluster_fe_error_mean"]
        elseif haskey(data, "cluster_pf_error_mean")
            return data["cluster_pf_error_mean"]
        elseif haskey(data, "cluster_error_mean")
            return data["cluster_error_mean"]
        end
    elseif error_type == "loop"
        if haskey(data, "loop_fe_error_mean")
            return data["loop_fe_error_mean"]
        elseif haskey(data, "loop_pf_error_mean")
            return data["loop_pf_error_mean"]
        elseif haskey(data, "loop_error_mean")
            return data["loop_error_mean"]
        end
    end
    return nothing
end

function plot_1_combined_error_vs_eta(cluster_results::Dict, loop_results::Dict)
    """
    Plot 1: BP free energy per site error and cluster/loop-corrected free energy per site errors vs η for different weights.
    Shows how both correction methods improve free energy per site accuracy across different η values.
    Separate plots for each N. Handles missing data gracefully.
    """
    # Get all available N values from both datasets
    all_N = unique(vcat(collect(keys(cluster_results)), collect(keys(loop_results))))
    
    for N in sort(all_N)
        p = plot(title="Combined Corrections: Free Energy per Site Error vs η (N=$N)", 
                 xlabel="η", ylabel="Free Energy per Site Error", 
                 legend=:topright, size=(1000, 700), dpi=300)
        
        # Colors for different methods and weights
        bp_color = :red
        cluster_colors = [:darkblue, :darkgreen, :indigo, :darkorange, :purple, :brown]
        loop_colors = [:blue, :green, :purple, :orange, :magenta, :cyan]
        markers = [:circle, :square, :diamond, :utriangle, :star5, :hexagon]
        
        # Collect all available η values from both datasets
        all_etas = []
        
        # Get η values from cluster data
        if haskey(cluster_results, N)
            for w_data in values(cluster_results[N])
                append!(all_etas, collect(keys(w_data)))
            end
        end
        
        # Get η values from loop data
        if haskey(loop_results, N)
            for w_data in values(loop_results[N])
                append!(all_etas, collect(keys(w_data)))
            end
        end
        
        all_etas = sort(unique(all_etas))
        
        if isempty(all_etas)
            println("No η values available for N=$N, skipping...")
            continue
        end
        
        # Plot BP error (same for both methods, use whichever is available)
        bp_means = []
        bp_etas = []
        
        for η in all_etas
            bp_err = nothing
            
            # Try to get BP error from cluster data first
            if haskey(cluster_results, N)
                for w in keys(cluster_results[N])
                    if haskey(cluster_results[N][w], η)
                        bp_err = get_error_value(cluster_results[N][w][η], "bp")
                        if bp_err !== nothing
                            break
                        end
                    end
                end
            end
            
            # If not found, try loop data
            if bp_err === nothing && haskey(loop_results, N)
                for w in keys(loop_results[N])
                    if haskey(loop_results[N][w], η)
                        bp_err = get_error_value(loop_results[N][w][η], "bp")
                        if bp_err !== nothing
                            break
                        end
                    end
                end
            end
            
            if bp_err !== nothing
                push!(bp_means, bp_err)
                push!(bp_etas, η)
            end
        end
        
        if !isempty(bp_means)
            plot!(p, bp_etas, bp_means, label="BP (no correction)", 
                  color=bp_color, linewidth=3, marker=:circle, markersize=5)
        end
        
        # Plot cluster corrections
        if haskey(cluster_results, N)
            weights = sort(collect(keys(cluster_results[N])))
            
            for (i, w) in enumerate(weights)
                # Collect data points that exist
                η_available = []
                cluster_means = []
                
                for η in all_etas
                    if haskey(cluster_results[N][w], η)
                        cluster_err = get_error_value(cluster_results[N][w][η], "cluster")
                        if cluster_err !== nothing
                            push!(η_available, η)
                            push!(cluster_means, cluster_err)
                        end
                    end
                end
                
                if !isempty(cluster_means)
                    color_idx = min(i, length(cluster_colors))
                    marker_idx = min(i, length(markers))
                    
                    # Plot cluster corrections (solid lines)
                    plot!(p, η_available, cluster_means, 
                          label="Cluster w=$w", color=cluster_colors[color_idx], linewidth=2, 
                          marker=markers[marker_idx], markersize=4, linestyle=:solid)
                end
            end
        end
        
        # Plot loop corrections
        if haskey(loop_results, N)
            weights = sort(collect(keys(loop_results[N])))
            
            for (i, w) in enumerate(weights)
                # Collect data points that exist
                η_available = []
                loop_means = []
                
                for η in all_etas
                    if haskey(loop_results[N][w], η)
                        loop_err = get_error_value(loop_results[N][w][η], "loop")
                        if loop_err !== nothing
                            push!(η_available, η)
                            push!(loop_means, loop_err)
                        end
                    end
                end
                
                if !isempty(loop_means)
                    color_idx = min(i, length(loop_colors))
                    marker_idx = min(i, length(markers))
                    
                    # Plot loop corrections (dotted lines)
                    plot!(p, η_available, loop_means, 
                          label="Loop w=$w", color=loop_colors[color_idx], linewidth=2, 
                          marker=markers[marker_idx], markersize=4, linestyle=:dot)
                end
            end
        end
        
        # Add log scale for y-axis
        plot!(p, yscale=:log10)
        
        savefig(p, "combined_fe_error_vs_eta_N$N.png")
        display(p)
        println("Saved: combined_fe_error_vs_eta_N$N.png")
    end
end

function plot_2_combined_error_vs_weight(cluster_results::Dict, loop_results::Dict, η_array::Vector{Float64}=[0.1, 0.5, 0.8])
    """
    Plot 2: Free energy per site error vs weight order for multiple η values.
    Shows both cluster and loop corrections for different η values on the same plot.
    Cluster corrections: solid lines
    Loop corrections: dotted lines
    """
    # Get all available N values from both datasets
    all_N = unique(vcat(collect(keys(cluster_results)), collect(keys(loop_results))))
    
    for N in sort(all_N)
        # Find all available η values across all weights from both datasets
        all_etas = []
        
        if haskey(cluster_results, N)
            for w_data in values(cluster_results[N])
                append!(all_etas, collect(keys(w_data)))
            end
        end
        
        if haskey(loop_results, N)
            for w_data in values(loop_results[N])
                append!(all_etas, collect(keys(w_data)))
            end
        end
        
        all_etas = sort(unique(all_etas))
        
        if isempty(all_etas)
            println("No η values available for N=$N, skipping...")
            continue
        end
        
        # Create plot for this N
        p = plot(title="Combined Corrections: Free Energy per Site Error vs Weight (N=$N)", 
                 xlabel="Weight", ylabel="Free Energy per Site Error",
                 yscale=:log,
                 legend=:topright, size=(1000, 700), dpi=300)
        
        # Colors and markers for different η values
        colors = [:red, :blue, :green, :purple, :orange, :brown, :pink, :gray]
        markers = [:circle, :square, :diamond, :utriangle, :star5, :hexagon, :cross, :plus]
        
        # Plot for each η value in the array
        for (i, η_target) in enumerate(η_array)
            # Find closest available η
            η_actual = all_etas[argmin(abs.(all_etas .- η_target))]
            
            color = colors[mod1(i, length(colors))]
            marker = markers[mod1(i, length(markers))]
            
            # Plot cluster corrections
            if haskey(cluster_results, N)
                weights_cluster = Float64[]
                errors_cluster = Float64[]
                
                # Get BP error
                bp_err = nothing
                for w in keys(cluster_results[N])
                    if haskey(cluster_results[N][w], η_actual)
                        bp_err = get_error_value(cluster_results[N][w][η_actual], "bp")
                        if bp_err !== nothing
                            break
                        end
                    end
                end
                
                if bp_err !== nothing
                    # Add w=0 (BP error) first
                    push!(weights_cluster, 0.0)
                    push!(errors_cluster, bp_err)
                    
                    # Add cluster corrections for other weights
                    for w in sort(collect(keys(cluster_results[N])))
                        if w == 0
                            continue  # Already handled BP error above
                        end
                        if haskey(cluster_results[N][w], η_actual)
                            cluster_err = get_error_value(cluster_results[N][w][η_actual], "cluster")
                            if cluster_err !== nothing
                                push!(weights_cluster, w)
                                push!(errors_cluster, cluster_err)
                            end
                        end
                    end
                    
                    if length(weights_cluster) > 1
                        # Plot cluster curve (solid line)
                        plot!(p, weights_cluster, errors_cluster, 
                              label="Cluster η=$η_actual", 
                              color=color, 
                              linewidth=2, 
                              marker=marker, 
                              markersize=6,
                              linestyle=:solid)
                    end
                end
            end
            
            # Plot loop corrections
            if haskey(loop_results, N)
                weights_loop = Float64[]
                errors_loop = Float64[]
                
                # Get BP error
                bp_err = nothing
                for w in keys(loop_results[N])
                    if haskey(loop_results[N][w], η_actual)
                        bp_err = get_error_value(loop_results[N][w][η_actual], "bp")
                        if bp_err !== nothing
                            break
                        end
                    end
                end
                
                if bp_err !== nothing
                    # Add w=0 (BP error) first
                    push!(weights_loop, 0.0)
                    push!(errors_loop, bp_err)
                    
                    # Add loop corrections for other weights
                    for w in sort(collect(keys(loop_results[N])))
                        if w == 0
                            continue  # Already handled BP error above
                        end
                        if haskey(loop_results[N][w], η_actual)
                            loop_err = get_error_value(loop_results[N][w][η_actual], "loop")
                            if loop_err !== nothing
                                push!(weights_loop, w)
                                push!(errors_loop, loop_err)
                            end
                        end
                    end
                    
                    if length(weights_loop) > 1
                        # Plot loop curve (dotted line)
                        plot!(p, weights_loop, errors_loop, 
                              label="Loop η=$η_actual", 
                              color=color, 
                              linewidth=2, 
                              marker=marker, 
                              markersize=6,
                              linestyle=:dot)
                    end
                end
            end
        end
        
        # Save the plot
        filename = "combined_error_vs_weight_multi_eta_N$(N).png"
        savefig(p, filename)
        display(p)
        println("Saved: $filename")
    end
end

# Main execution
function main()
    println("="^60)
    println("🔄 Combined Cluster & Loop Correction Plotting Script")
    println("="^60)
    
    println("Loading cluster correction results...")
    cluster_results = load_cluster_results()
    
    println("Loading loop correction results...")
    loop_results = load_loop_results()
    
    if isempty(cluster_results) && isempty(loop_results)
        println("❌ No correction results found!")
        println("Make sure cluster_correction_results and/or loop_correction_results directories exist with .jld2 files.")
        return
    end
    
    println("\n📊 Available data:")
    
    if !isempty(cluster_results)
        println("\nCluster correction data:")
        for N in sort(collect(keys(cluster_results)))
            weights = sort(collect(keys(cluster_results[N])))
            total_combinations = sum([length(keys(cluster_results[N][w])) for w in weights])
            println("  N=$N: weights=$weights, total combinations=$total_combinations")
        end
    end
    
    if !isempty(loop_results)
        println("\nLoop correction data:")
        for N in sort(collect(keys(loop_results)))
            weights = sort(collect(keys(loop_results[N])))
            total_combinations = sum([length(keys(loop_results[N][w])) for w in weights])
            println("  N=$N: weights=$weights, total combinations=$total_combinations")
        end
    end
    
    println("\n🎨 Generating combined plots...")
    
    println("\nGenerating Plot 1: Combined Free Energy per Site Error vs η...")
    plot_1_combined_error_vs_eta(cluster_results, loop_results)
    
    println("\nGenerating Plot 2: Combined Free Energy per Site Error vs weight...")
    plot_2_combined_error_vs_weight(cluster_results, loop_results, [0.1, 0.3, 0.5, 0.7, 0.9])
    
    println("\n✅ Done! Check the generated PNG files.")
    println("Files generated:")
    println("  - combined_fe_error_vs_eta_N*.png")
    println("  - combined_error_vs_weight_multi_eta_N*.png")
    println("\nLegend:")
    println("  🔹 Cluster corrections: solid lines")
    println("  🔸 Loop corrections: dotted lines")
    println("  🔴 BP (no correction): red solid line")
end

# Run the script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
