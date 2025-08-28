using JLD2, Plots, Statistics
using Dates

"""
Collect and organize results from the Gallager BP SLURM job array.
This script should be run after all SLURM jobs have completed.
"""

function collect_gallager_results(results_dir="gallagerBP0results")
    # Get all result files
    files = filter(x -> startswith(x, "results_task") && endswith(x, ".jld2"), readdir(results_dir))
    
    if isempty(files)
        @warn "No result files found in $results_dir"
        return
    end
    
    # Extract parameters from first file to get structure
    first_file = joinpath(results_dir, files[1])
    data = load(first_file)
    
    # Get unique n values and p values
    n_values = Set{Int}()
    p_values = Set{Float64}()
    
    for file in files
        filepath = joinpath(results_dir, file)
        try
            data = load(filepath)
            push!(n_values, data["n"])
            push!(p_values, data["p"])
        catch e
            @warn "Could not load $file: $e"
        end
    end
    
    n_values = sort(collect(n_values))
    p_values = sort(collect(p_values))
    
    @info "Found $(length(files)) result files"
    @info "n values: $n_values"
    @info "p values: $(length(p_values)) points from $(minimum(p_values)) to $(maximum(p_values))"
    
    # Initialize result arrays
    results = Dict()
    for n in n_values
        results[n] = Dict(
            "p_values" => Float64[],
            "mean_error_rate_loops" => Float64[],
            "std_error_rate_loops" => Float64[],
            "mean_error_rate_no_loops" => Float64[],
            "std_error_rate_no_loops" => Float64[],
        )
    end
    
    # Collect all results
    for file in files
        filepath = joinpath(results_dir, file)
        try
            data = load(filepath)
            n = data["n"]
            
            push!(results[n]["p_values"], data["p"])
            push!(results[n]["mean_error_rate_loops"], data["mean_error_rate_loops"])
            push!(results[n]["std_error_rate_loops"], data["std_error_rate_loops"])
            push!(results[n]["mean_error_rate_no_loops"], data["mean_error_rate_no_loops"])
            push!(results[n]["std_error_rate_no_loops"], data["std_error_rate_no_loops"])
        catch e
            @warn "Could not process $file: $e"
        end
    end
    
    # Sort results by p value for each n
    for n in n_values
        perm = sortperm(results[n]["p_values"])
        for key in keys(results[n])
            results[n][key] = results[n][key][perm]
        end
    end
    
    # Save collected results with timestamp
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    collected_filename = joinpath(results_dir, "collected_results_$(timestamp).jld2")
    
    # Get metadata from first file
    first_data = load(joinpath(results_dir, files[1]))
    metadata = Dict(
        "d_v" => first_data["d_v"],
        "d_c" => first_data["d_c"],
        "num_batches" => first_data["num_batches"],
        "samples_per_batch" => first_data["samples_per_batch"],
        "max_loop_order" => first_data["max_loop_order"],
        "timestamp" => timestamp,
        "num_files_processed" => length(files)
    )
    
    save(collected_filename, 
         "results", results,
         "n_values", n_values,
         "metadata", metadata)
    
    @info "Saved collected results to $collected_filename"
    
    # Create summary plots
    create_summary_plots(results, n_values, timestamp, results_dir)
    
    return results, n_values, metadata
end

function create_summary_plots(results, n_values, timestamp, results_dir)
    @info "Creating summary plots..."
    
    # Plot error rates vs p for different n values
    p1 = plot(title="Error Rate vs Physical Error Rate p", 
              xlabel="Physical Error Rate p", 
              ylabel="Logical Error Rate",
              yscale=:log10,
              legend=:bottomright)
    
    colors = [:blue, :red, :green, :orange, :purple]
    
    for (i, n) in enumerate(n_values)
        color = colors[mod1(i, length(colors))]
        
        # Plot with loops
        plot!(p1, results[n]["p_values"], results[n]["mean_error_rate_loops"],
              ribbon=results[n]["std_error_rate_loops"],
              label="n=$n (with loops)", color=color, linestyle=:solid,
              marker=:circle, markersize=3)
        
        # Plot without loops
        plot!(p1, results[n]["p_values"], results[n]["mean_error_rate_no_loops"],
              ribbon=results[n]["std_error_rate_no_loops"],
              label="n=$n (no loops)", color=color, linestyle=:dash,
              marker=:square, markersize=3)
    end
    
    savefig(p1, joinpath(results_dir, "error_rates_vs_p_$(timestamp).png"))
    savefig(p1, joinpath(results_dir, "error_rates_vs_p_$(timestamp).pdf"))
    
    # Plot improvement due to loops
    p2 = plot(title="Loop Correction Improvement", 
              xlabel="Physical Error Rate p", 
              ylabel="Error Rate Ratio (no loops / with loops)",
              yscale=:log10,
              legend=:topright)
    
    for (i, n) in enumerate(n_values)
        color = colors[mod1(i, length(colors))]
        ratio = results[n]["mean_error_rate_no_loops"] ./ results[n]["mean_error_rate_loops"]
        plot!(p2, results[n]["p_values"], ratio,
              label="n=$n", color=color, marker=:circle, markersize=3)
    end
    
    hline!(p2, [1.0], color=:black, linestyle=:dash, label="No improvement")
    
    savefig(p2, joinpath(results_dir, "loop_improvement_$(timestamp).png"))
    savefig(p2, joinpath(results_dir, "loop_improvement_$(timestamp).pdf"))
    
    @info "Plots saved with timestamp: $timestamp"
end

# Run collection if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    results, n_values, metadata = collect_gallager_results()
    
    @info "Collection complete!"
    @info "Results structure:"
    for n in n_values
        @info "  n=$n: $(length(results[n]["p_values"])) data points"
    end
end
