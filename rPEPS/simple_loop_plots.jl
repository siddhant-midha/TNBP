using JLD2, FileIO, Plots, Statistics
using Printf

"""
Simple plotting script for loop correction analysis.
Generates free energy per site error comparison plots:
1. BP free energy per site error vs η compared with loop-corrected errors for different weights
2. Free energy per site error vs weight order for fixed η
3. Improvement summary across different system sizes
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

function load_loop_results(results_dir::String)
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

function plot_1_loop_error_vs_eta(loop_results::Dict)
    mode = get_mode()
    suffix = mode
    figs_dir = ensure_figs_dir()
    for N in sort(collect(keys(loop_results)))
        p = plot(title="Loop Correction: Free Energy per Site Error vs η (N=$N)", 
                 xlabel="η", ylabel="Free Energy per Site Error", 
                 legend=:topright, size=(900, 600), dpi=300)
        
        # Get available weights for this N
        weights = sort(collect(keys(loop_results[N])))
        
        if isempty(weights)
            println("No weights available for N=$N, skipping...")
            continue
        end
        
        # Get η values from the first available weight
        sample_w = weights[1]
        if isempty(loop_results[N][sample_w])
            println("No η values available for N=$N, w=$sample_w, skipping...")
            continue
        end
        
        η_vals = sort(collect(keys(loop_results[N][sample_w])))
        
        # Colors and markers for different correction types
        bp_color = :red
        loop_colors = [:darkblue, :darkgreen, :indigo, :darkorange, :purple]
        markers = [:circle, :square, :diamond, :utriangle, :star5]
        
        # Plot BP error (same for all weights, use first available)
        bp_means = []
        for η in η_vals
            if haskey(loop_results[N][sample_w], η)
                # Use new free energy per site error keys with fallback to old keys
                data = loop_results[N][sample_w][η]
                if haskey(data, "bp_fe_error_mean")
                    bp_error = data["bp_fe_error_mean"]
                elseif haskey(data, "bp_error_mean") 
                    bp_error = data["bp_error_mean"]
                else
                    println("Warning: No BP error data found for N=$N, w=$sample_w, η=$η")
                    continue
                end
                push!(bp_means, bp_error)
            end
        end
        
        if !isempty(bp_means)
            plot!(p, η_vals[1:length(bp_means)], bp_means, label="BP (no correction)", 
                  color=bp_color, linewidth=2, marker=:circle, markersize=4)
        end
        
        # Plot loop corrections for each available weight
        for (i, w) in enumerate(weights)
            # Collect data points that exist
            η_available = []
            loop_means = []
            
            for η in η_vals
                if haskey(loop_results[N][w], η)
                    push!(η_available, η)
                    # Use new free energy per site error keys with fallback to old keys
                    data = loop_results[N][w][η]
                    if haskey(data, "loop_fe_error_mean")
                        loop_error = data["loop_fe_error_mean"]
                    elseif haskey(data, "loop_error_mean")
                        loop_error = data["loop_error_mean"]
                    else
                        println("Warning: No loop error data found for N=$N, w=$w, η=$η")
                        continue
                    end
                    push!(loop_means, loop_error)
                end
            end
            
            if !isempty(loop_means)
                color_idx = min(i, length(loop_colors))
                marker_idx = min(i, length(markers))
                
                # Plot loop corrections
                plot!(p, η_available, loop_means, 
                      label="Loop w=$w", color=loop_colors[color_idx], linewidth=2, 
                      marker=markers[marker_idx], markersize=4, linestyle=:dash)
            end
        end
        
        # Add log scale for y-axis
        plot!(p, yscale=:log10)
        
        filename = joinpath(figs_dir, "loop_fe_error_vs_eta_N$(N)_$(suffix).png")
        savefig(p, filename)
        display(p)
        println("Saved: $(filename)")
    end
end

function plot_2_loop_error_vs_weight(loop_results::Dict, η_array::Vector{Float64}=[0.1, 0.5, 0.8])
    mode = get_mode()
    suffix = mode
    figs_dir = ensure_figs_dir()
    for N in sort(collect(keys(loop_results)))
        # Find all available η values across all weights
        all_etas = []
        for w_data in values(loop_results[N])
            append!(all_etas, collect(keys(w_data)))
        end
        all_etas = unique(all_etas)
        
        if isempty(all_etas)
            println("No η values available for N=$N, skipping...")
            continue
        end
        
        # Create plot for this N
        p = plot(title="Loop Correction: Free Energy per Site Error vs Weight (N=$N)", 
                 xlabel="Cluster Weight", ylabel="Free Energy per Site Error",
                 yscale=:log,
                 legend=:topright, size=(900, 600), dpi=300)
        
        # Colors and markers for different η values
        colors = [:red, :blue, :green, :purple, :orange, :brown, :pink, :gray]
        markers = [:circle, :square, :diamond, :utriangle, :star5, :hexagon, :cross, :plus]
        
        # Plot for each η value in the array
        for (i, η_target) in enumerate(η_array)
            # Find closest available η
            η_actual = all_etas[argmin(abs.(all_etas .- η_target))]
            
            # Collect data points for available weights - single curve per η
            weights = Float64[]
            errors = Float64[]
            
            # Get BP error from any available weight that has this η
            bp_err = nothing
            for w in keys(loop_results[N])
                if haskey(loop_results[N][w], η_actual)
                    # Use new free energy per site error keys with fallback to old keys
                    data = loop_results[N][w][η_actual]
                    if haskey(data, "bp_fe_error_mean")
                        bp_err = data["bp_fe_error_mean"]
                    elseif haskey(data, "bp_error_mean")
                        bp_err = data["bp_error_mean"]
                    end
                    if bp_err !== nothing
                        break
                    end
                end
            end
            
            if bp_err === nothing
                println("No BP error data available for N=$N, η=$η_actual, skipping...")
                continue
            end
            
            # Add w=0 (BP error) first
            push!(weights, 0.0)  # Use actual weight 0 for BP
            push!(errors, bp_err)
            
            # Add loop corrections for other weights
            for w in sort(collect(keys(loop_results[N])))
                if w == 0
                    continue  # Already handled BP error above
                end
                if haskey(loop_results[N][w], η_actual)
                    # Use new free energy per site error keys with fallback to old keys
                    data = loop_results[N][w][η_actual]
                    if haskey(data, "loop_fe_error_mean")
                        loop_error = data["loop_fe_error_mean"]
                    elseif haskey(data, "loop_error_mean")
                        loop_error = data["loop_error_mean"]
                    else
                        println("Warning: No loop error data found for N=$N, w=$w, η=$η_actual")
                        continue
                    end
                    push!(weights, w)
                    push!(errors, loop_error)
                end
            end
            
            if length(weights) <= 1
                println("Insufficient data for N=$N, η=$η_actual")
                continue
            end
            
            # Plot single curve for this η value
            color = colors[mod1(i, length(colors))]
            marker = markers[mod1(i, length(markers))]
            
            # Plot single curve showing error vs weight (BP at w=0.1, loop corrections at higher w)
            plot!(p, weights, errors, 
                  label="η=$η_actual", 
                  color=color, 
                  linewidth=2, 
                  marker=marker, 
                  markersize=6)
        end
        
        # Save the plot
        filename = joinpath(figs_dir, "loop_error_vs_weight_multi_eta_N$(N)_$(suffix).png")
        savefig(p, filename)
        println("Saved: $(filename)")
    end
end

function plot_improvement_summary(loop_results::Dict)
    mode = get_mode()
    suffix = mode
    figs_dir = ensure_figs_dir()
    for N in sort(collect(keys(loop_results)))
        p = plot(title="Loop Correction Improvement (N=$N)", 
                 xlabel="η", ylabel="Free Energy per Site Error Improvement (BP - Loop)",
                 legend=:topright, size=(900, 600), dpi=300)
        
        # Collect available weights
        weights = sort(collect(keys(loop_results[N])))
        colors = [:darkblue, :darkgreen, :indigo, :darkorange, :purple]
        
        for (i, w) in enumerate(weights)
            η_vals = sort(collect(keys(loop_results[N][w])))
            
            loop_improvements = []
            η_available = []
            
            for η in η_vals
                loop_data = loop_results[N][w][η]
                # Use new free energy per site error keys with fallback to old keys
                if haskey(loop_data, "bp_fe_error_mean") && haskey(loop_data, "loop_fe_error_mean")
                    bp_err = loop_data["bp_fe_error_mean"]
                    loop_err = loop_data["loop_fe_error_mean"]
                elseif haskey(loop_data, "bp_error_mean") && haskey(loop_data, "loop_error_mean")
                    bp_err = loop_data["bp_error_mean"]
                    loop_err = loop_data["loop_error_mean"]
                else
                    println("Warning: Missing error data for N=$N, w=$w, η=$η")
                    continue
                end
                loop_improvement = bp_err - loop_err
                
                push!(loop_improvements, loop_improvement)
                push!(η_available, η)
            end
            
            # Plot loop improvements
            if !isempty(loop_improvements)
                color_idx = min(i, length(colors))
                plot!(p, η_available, loop_improvements, 
                      label="Loop w=$w", linewidth=2, marker=:circle, markersize=3,
                      color=colors[color_idx])
            end
        end
        
        # Add horizontal line at y=0 for reference
        plot!(p, [minimum([minimum(collect(keys(loop_results[N][w]))) for w in weights]),
                   maximum([maximum(collect(keys(loop_results[N][w]))) for w in weights])], 
              [0, 0], color=:black, linestyle=:dot, alpha=0.5, label="No improvement")
        
        filename = joinpath(figs_dir, "loop_fe_improvement_summary_N$(N)_$(suffix).png")
        savefig(p, filename)
        display(p)
        println("Saved: $(filename)")
    end
end

# Main execution
function main()
    mode = get_mode()
    results_dir = mode == "complex" ? "loop_correction_results_complex" : "loop_correction_results"
    println("Loading loop correction results from $results_dir...")
    loop_results = load_loop_results(results_dir)
    ensure_figs_dir()
    
    if isempty(loop_results)
        println("❌ No loop correction results found!")
        println("Make sure loop_correction_results directory exists with .jld2 files.")
        return
    end
    
    println("\n📊 Available loop correction data:")
    for N in sort(collect(keys(loop_results)))
        weights = sort(collect(keys(loop_results[N])))
        total_combinations = sum([length(keys(loop_results[N][w])) for w in weights])
        println("  N=$N: weights=$weights, total combinations=$total_combinations")
        
        # Show success/failure pattern
        for w in weights
            η_vals = sort(collect(keys(loop_results[N][w])))
            println("    w=$w: η values = $η_vals")
        end
    end
    
    println("\n🎨 Generating plots...")
    
    println("\nGenerating Plot 1: Loop Free Energy per Site Error vs η...")
    plot_1_loop_error_vs_eta(loop_results)
    
    println("\nGenerating Plot 2: Loop Free Energy per Site Error vs weight...")
    plot_2_loop_error_vs_weight(loop_results, [0.1, 0.3, 0.5,0.7, 0.9, 1.0])
    
    println("\nGenerating Plot 3: Loop improvement summary...")
    plot_improvement_summary(loop_results)
    
    println("\n✅ Done! Check the generated PNG files.")
    println("Files generated:")
    println("  - loop_fe_error_vs_eta_N*_complex.png")
    println("  - loop_error_vs_weight_multi_eta_N*_complex.png") 
    println("  - loop_fe_improvement_summary_N*_complex.png")
end

# Run the script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
