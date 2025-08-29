using JLD2, FileIO, Plots, Statistics
using Printf

"""
Simple plotting script for cluster correction analysis.
Generates free energy per site error comparison plots:
1. BP free energy per site error vs η compared with cluster-corrected errors for different weights
2. Free energy per site error vs weight order for fixed η
Compatible with old partition function errors and new free energy per site errors.
"""

# Parse argument for mode
function get_mode()
    mode = "positive" # default
    for arg in ARGS
        if arg == "--complex"
            mode = "complex"
        elseif arg == "--positive"
            mode = "positive"
        end
    end
    return mode
end

function load_results_simple(results_dir::String="cluster_correction_results_complex")
    """
    Load and organize results from JLD2 files.
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
            println("Warning: Failed to load $file: $e")
        end
    end
    
    return results
end

function ensure_figs_dir()
    figs_dir = "figs"
    if !isdir(figs_dir)
        mkpath(figs_dir)
    end
    return figs_dir
end

function plot_1_error_vs_eta(results::Dict)
    """
    Plot 1: BP free energy per site error and cluster-corrected free energy per site errors vs η for different weights.
    Shows how cluster correction improves free energy per site accuracy across different η values.
    Separate plots for each N. Handles missing data gracefully.
    """
    mode = get_mode()
    suffix = mode
    figs_dir = ensure_figs_dir()
    for N in sort(collect(keys(results)))
        p = plot(title="Free Energy per Site Error vs η (N=$N)", xlabel="η", ylabel="Free Energy per Site Error", 
                 legend=:topright, size=(800, 600), dpi=300)
        
        # Get available weights for this N
        weights = sort(collect(keys(results[N])))
        
        if isempty(weights)
            println("No weights available for N=$N, skipping...")
            continue
        end
        
        # Get all unique η values across all weights (some may be missing for certain weights)
        all_etas = Set{Float64}()
        for w in weights
            union!(all_etas, collect(keys(results[N][w])))
        end
        η_vals = sort(collect(all_etas))
        
        if isempty(η_vals)
            println("No η values available for N=$N, skipping...")
            continue
        end
        
        # Plot BP error (use first available weight for BP data)
        sample_w = weights[1]
        bp_means = []
        η_bp_available = []
        
        for η in η_vals
            if haskey(results[N][sample_w], η)
                # Use new free energy per site error keys with fallback to old keys
                data = results[N][sample_w][η]
                if haskey(data, "bp_fe_error_mean")
                    bp_error = data["bp_fe_error_mean"]
                elseif haskey(data, "bp_pf_error_mean")
                    bp_error = data["bp_pf_error_mean"]
                elseif haskey(data, "bp_error_mean")
                    bp_error = data["bp_error_mean"]
                else
                    println("Warning: No BP error data found for N=$N, w=$sample_w, η=$η")
                    continue
                end
                push!(bp_means, bp_error)
                push!(η_bp_available, η)
            end
        end
        
        if !isempty(bp_means)
            plot!(p, η_bp_available, bp_means, label="BP (no correction)", 
                  color=:red, linewidth=2, marker=:circle, markersize=4)
        end
        
        # Plot cluster corrections for each weight
        colors = [:blue, :green, :purple, :orange, :brown]
        markers = [:square, :diamond, :utriangle, :star5, :hexagon]
        
        for (i, w) in enumerate(weights)
            cluster_means = []
            η_cluster_available = []
            
            for η in η_vals
                if haskey(results[N][w], η)
                    # Use new free energy per site error keys with fallback to old keys
                    data = results[N][w][η]
                    if haskey(data, "cluster_fe_error_mean")
                        cluster_error = data["cluster_fe_error_mean"]
                    elseif haskey(data, "cluster_pf_error_mean")
                        cluster_error = data["cluster_pf_error_mean"]
                    elseif haskey(data, "cluster_error_mean")
                        cluster_error = data["cluster_error_mean"]
                    else
                        println("Warning: No cluster error data found for N=$N, w=$w, η=$η")
                        continue
                    end
                    push!(cluster_means, cluster_error)
                    push!(η_cluster_available, η)
                end
            end
            
            if !isempty(cluster_means)
                color = i <= length(colors) ? colors[i] : :black
                marker = i <= length(markers) ? markers[i] : :circle
                plot!(p, η_cluster_available, cluster_means, 
                      label="Cluster w=$w", color=color, linewidth=2, 
                      marker=marker, markersize=4)
            else
                println("Warning: No data available for N=$N, w=$w")
            end
        end
        
        # Add log scale for y-axis
        plot!(p, yscale=:log10)
        
        filename = joinpath(figs_dir, "cluster_fe_error_vs_eta_N$(N)_$(suffix).png")
        savefig(p, filename)
        display(p)
        println("Saved: $(filename)")
    end
end

function plot_2_error_vs_weight(results::Dict, η_array::Vector{Float64}=[0.1, 0.5, 0.8])
    """
    Plot cluster-corrected errors as a function of cluster weight.
    Plots multiple η values on the same plot. Handles missing data gracefully.
    """
    mode = get_mode()
    suffix = mode
    figs_dir = ensure_figs_dir()
    for N in sort(collect(keys(results)))
        # Find all available η values across all weights
        all_etas = []
        for w_data in values(results[N])
            append!(all_etas, collect(keys(w_data)))
        end
        all_etas = unique(all_etas)
        
        if isempty(all_etas)
            println("No η values available for N=$N, skipping...")
            continue
        end
        
        # Create plot for this N
        p = plot(title="Free Energy per Site Error vs Cluster Weight (N=$N)", 
                 xlabel="Cluster Weight", ylabel="Free Energy per Site Error",
                 yscale=:log,
                 legend=:topright, size=(800, 600), dpi=300)
        
        # Colors and markers for different η values
        colors = [:red, :blue, :green, :purple, :orange, :brown, :pink, :gray]
        markers = [:circle, :square, :diamond, :utriangle, :star5, :hexagon, :cross, :plus]
        
        # Plot for each η value in the array
        for (i, η_target) in enumerate(η_array)
            # Find closest available η
            η_actual = all_etas[argmin(abs.(all_etas .- η_target))]
            
            # Collect data points for available weights
            weights = Float64[]
            cluster_errors = Float64[]
            
            # Get BP error (weight 0 first)
            bp_err = nothing
            if haskey(results[N], 0) && haskey(results[N][0], η_actual)
                data = results[N][0][η_actual]
                if haskey(data, "bp_fe_error_mean")
                    bp_err = data["bp_fe_error_mean"]
                elseif haskey(data, "bp_pf_error_mean")
                    bp_err = data["bp_pf_error_mean"]
                elseif haskey(data, "bp_error_mean")
                    bp_err = data["bp_error_mean"]
                end
            end
            
            # If no weight 0 data, try to get BP error from any weight
            if bp_err === nothing
                for w in sort(collect(keys(results[N])))
                    if haskey(results[N][w], η_actual)
                        data = results[N][w][η_actual]
                        if haskey(data, "bp_fe_error_mean")
                            bp_err = data["bp_fe_error_mean"]
                        elseif haskey(data, "bp_pf_error_mean")
                            bp_err = data["bp_pf_error_mean"]
                        elseif haskey(data, "bp_error_mean")
                            bp_err = data["bp_error_mean"]
                        end
                        if bp_err !== nothing
                            break
                        end
                    end
                end
            end
            
            if bp_err === nothing
                println("No BP error data available for N=$N, η=$η_actual, skipping...")
                continue
            end
            
            # Add BP error at weight 0 (actual weight)
            push!(weights, 0.0)
            push!(cluster_errors, bp_err)
        
            # Add cluster weights that have data for this η
            for w in sort(collect(keys(results[N])))
                if w == 0
                    continue  # Already handled BP error above
                end
                if haskey(results[N][w], η_actual)
                    # Use new free energy per site error keys with fallback to old keys
                    data = results[N][w][η_actual]
                    if haskey(data, "cluster_fe_error_mean")
                        cl_err = data["cluster_fe_error_mean"]
                    elseif haskey(data, "cluster_pf_error_mean")
                        cl_err = data["cluster_pf_error_mean"]
                    elseif haskey(data, "cluster_error_mean")
                        cl_err = data["cluster_error_mean"]
                    else
                        println("Warning: No cluster error data found for N=$N, w=$w, η=$η_actual")
                        continue
                    end
                    push!(weights, w)
                    push!(cluster_errors, cl_err)
                end
            end
            
            if length(weights) <= 1
                println("Insufficient data for N=$N, η=$η_actual")
                continue
            end
            
            # Plot this η series
            color = colors[mod1(i, length(colors))]
            marker = markers[mod1(i, length(markers))]
            
            plot!(p, weights, cluster_errors, 
                  label="η=$η_actual", 
                  color=color, 
                  marker=marker, 
                  markersize=6, 
                  linewidth=2)
        end
        
        # Save the plot
        filename = joinpath(figs_dir, "cluster_error_vs_weight_multi_eta_N$(N)_$(suffix).png")
        savefig(p, filename)
        println("Saved plot: $(filename)")
    end
end

# Main execution
function main()
    mode = get_mode()
    results_dir = mode == "complex" ? "cluster_correction_results_complex" : "cluster_correction_results_positive"
    println("Loading results from $results_dir...")
    results = load_results_simple(results_dir)
    ensure_figs_dir()
    
    if isempty(results)
        println("❌ No cluster correction results found!")
        println("Make sure cluster_correction_results directory exists with .jld2 files.")
        return
    end
    
    println("\n📊 Available data:")
    for N in sort(collect(keys(results)))
        weights = sort(collect(keys(results[N])))
        total_combinations = sum([length(keys(results[N][w])) for w in weights])
        println("  N=$N: weights=$weights, total combinations=$total_combinations")
        
        # Show data availability pattern
        for w in weights
            η_vals = sort(collect(keys(results[N][w])))
            println("    w=$w: η values = $η_vals")
        end
    end
    
    println("\n🎨 Generating plots...")
    
    println("\nGenerating Plot 1: Cluster Free Energy per Site Error vs η...")
    plot_1_error_vs_eta(results)
    
    println("\nGenerating Plot 2: Cluster Free Energy per Site Error vs weight...")
    plot_2_error_vs_weight(results, [0.1,0.2,0.3])  # Multiple η values
    
    println("\n✅ Done! Check the generated PNG files.")
    println("Files generated:")
    println("  - cluster_fe_error_vs_eta_N*_complex.png")
    println("  - cluster_error_vs_weight_multi_eta_N*_complex.png")
end

# Run the script
main()
