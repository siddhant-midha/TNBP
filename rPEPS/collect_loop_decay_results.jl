using JLD2, FileIO, Statistics
using Printf

"""
Script to collect and summarize loop decay results from SLURM job array runs.
"""

function collect_loop_decay_results(results_dir::String="loop_decay_results")
    """
    Collect all loop decay results and create a summary.
    """
    if !isdir(results_dir)
        println("❌ Directory $results_dir not found!")
        return
    end
    
    files = filter(f -> endswith(f, ".jld2"), readdir(results_dir))
    
    if isempty(files)
        println("❌ No .jld2 files found in $results_dir")
        return
    end
    
    println("📊 Found $(length(files)) result files")
    
    # Organize results by N and η
    organized_results = Dict()
    successful_runs = 0
    failed_runs = 0
    
    for file in files
        try
            filepath = joinpath(results_dir, file)
            data = load(filepath, "results")
            
            N = data["N"]
            η = data["eta"]
            nsamples = data["nsamples"]
            total_loops = data["total_loops_processed"]
            
            if !haskey(organized_results, N)
                organized_results[N] = Dict()
            end
            
            organized_results[N][η] = Dict(
                "nsamples" => nsamples,
                "total_loops_processed" => total_loops,
                "avg_contributions" => data["avg_contributions_by_weight"],
                "filename" => file
            )
            
            successful_runs += 1
            
        catch e
            println("⚠️  Failed to load $file: $e")
            failed_runs += 1
        end
    end
    
    println("\n📈 Summary:")
    println("  ✅ Successful runs: $successful_runs")
    println("  ❌ Failed runs: $failed_runs")
    
    # Print organized summary
    for N in sort(collect(keys(organized_results)))
        println("\n🔹 N = $N:")
        η_values = sort(collect(keys(organized_results[N])))
        
        for η in η_values
            data = organized_results[N][η]
            avg_contribs = data["avg_contributions"]
            weights = sort(collect(keys(avg_contribs)))
            
            println("  η = $η: $(data["nsamples"]) samples, $(data["total_loops_processed"]) total loops")
            println("    Weights with data: $weights")
            
            # Show some statistics
            if !isempty(weights)
                min_weight, max_weight = minimum(weights), maximum(weights)
                min_contrib = avg_contribs[min_weight]
                max_contrib = avg_contribs[max_weight]
                
                if min_contrib > 0 && max_contrib > 0
                    decay_factor = max_contrib / min_contrib
                    println("    Contribution range: $(max_contrib:.2e) (w=$max_weight) → $(min_contrib:.2e) (w=$min_weight)")
                    println("    Decay factor: $(decay_factor:.2e)")
                end
            end
        end
    end
    
    # Save organized results
    summary_file = joinpath(results_dir, "loop_decay_summary.jld2")
    save(summary_file, "organized_results", organized_results)
    println("\n💾 Saved organized results to: $summary_file")
    
    return organized_results
end

function check_missing_combinations()
    """
    Check which (N, η) combinations are missing from the results.
    """
    # Expected combinations
    N_values = [6, 8, 10]
    η_values = 0.1:0.1:2.0
    
    results = collect_loop_decay_results()
    
    println("\n🔍 Checking for missing combinations...")
    
    missing_count = 0
    total_expected = length(N_values) * length(η_values)
    
    for N in N_values
        for η in η_values
            if !haskey(results, N) || !haskey(results[N], η)
                println("❌ Missing: N=$N, η=$η")
                missing_count += 1
            end
        end
    end
    
    completed_count = total_expected - missing_count
    completion_rate = (completed_count / total_expected) * 100
    
    println("\n📊 Completion Status:")
    println("  Total expected: $total_expected")
    println("  Completed: $completed_count")
    println("  Missing: $missing_count")
    println("  Completion rate: $(round(completion_rate, digits=1))%")
    
    if missing_count > 0
        println("\n💡 To rerun missing combinations, you may need to:")
        println("  1. Check SLURM logs for failed jobs")
        println("  2. Resubmit specific array indices")
        println("  3. Or rerun the entire job array")
    end
end

function analyze_decay_patterns(results_dir::String="loop_decay_results")
    """
    Analyze patterns in the loop decay data.
    """
    println("\n🔬 Analyzing decay patterns...")
    
    results = collect_loop_decay_results(results_dir)
    
    if isempty(results)
        return
    end
    
    # Analyze decay rates
    decay_analysis = Dict()
    
    for N in sort(collect(keys(results)))
        decay_analysis[N] = Dict()
        
        for η in sort(collect(keys(results[N])))
            avg_contribs = results[N][η]["avg_contributions"]
            
            if length(avg_contribs) >= 3
                weights = sort(collect(keys(avg_contribs)))
                contributions = [avg_contribs[w] for w in weights]
                
                # Calculate exponential decay rate
                if all(contributions .> 0)
                    log_contributions = log.(contributions)
                    
                    # Linear fit: log(y) = a - b*x
                    n = length(weights)
                    sum_w = sum(weights)
                    sum_logc = sum(log_contributions)
                    sum_w2 = sum(weights.^2)
                    sum_w_logc = sum(weights .* log_contributions)
                    
                    decay_rate = -(n * sum_w_logc - sum_w * sum_logc) / (n * sum_w2 - sum_w^2)
                    intercept = (sum_logc + decay_rate * sum_w) / n
                    
                    # R-squared for fit quality
                    y_mean = mean(log_contributions)
                    ss_tot = sum((log_contributions .- y_mean).^2)
                    y_pred = intercept .- decay_rate .* weights
                    ss_res = sum((log_contributions .- y_pred).^2)
                    r_squared = 1 - ss_res / ss_tot
                    
                    decay_analysis[N][η] = Dict(
                        "decay_rate" => decay_rate,
                        "intercept" => intercept,
                        "r_squared" => r_squared,
                        "num_points" => n
                    )
                end
            end
        end
    end
    
    # Print decay analysis
    println("\n📉 Exponential Decay Analysis:")
    for N in sort(collect(keys(decay_analysis)))
        println("\nN = $N:")
        for η in sort(collect(keys(decay_analysis[N])))
            analysis = decay_analysis[N][η]
            println(@sprintf("  η = %.1f: decay_rate = %.3f, R² = %.3f (%d points)", 
                    η, analysis["decay_rate"], analysis["r_squared"], analysis["num_points"]))
        end
    end
    
    # Save decay analysis
    analysis_file = joinpath(results_dir, "decay_analysis.jld2")
    save(analysis_file, "decay_analysis", decay_analysis)
    println("\n💾 Saved decay analysis to: $analysis_file")
end

function main()
    println("="^60)
    println("🔄 Loop Decay Results Collection Script")
    println("="^60)
    
    # Collect and organize results
    collect_loop_decay_results()
    
    # Check for missing combinations
    check_missing_combinations()
    
    # Analyze decay patterns
    analyze_decay_patterns()
    
    println("\n✅ Collection and analysis complete!")
    println("\nNext steps:")
    println("  1. Run: julia plot_loop_decay.jl")
    println("  2. Check the generated plots")
    println("  3. Analyze decay patterns across different N and η values")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
