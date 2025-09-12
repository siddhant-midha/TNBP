# filepath: /Users/siddhantm/Desktop/TNBP/TNBP/aklt/aklt_ctmrg_benchmark.jl

using Plots
using Printf
include("aklt.jl")

"""
Benchmark CTMRG convergence for AKLT model by testing bond dimension and iteration effects.
"""

function benchmark_ctmrg_bond_dimension(;
                                       a1::Float64=sqrt(3/2), 
                                       a2::Float64=sqrt(6),
                                       χ_values=[50,100,150],
                                       nsteps::Int=500,
                                       cutoff::Float64=1e-18,
                                       save_plots::Bool=true,
                                       plot_dir::String="plots")
    """
    Benchmark CTMRG convergence vs bond dimension χ for fixed iteration number.
    
    Args:
        a1, a2: AKLT model parameters (default: sqrt(3/2), sqrt(6))
        χ_values: Bond dimensions to test
        nsteps: Fixed number of CTMRG iterations
        cutoff: Convergence cutoff
        save_plots: Whether to save plots
        plot_dir: Directory for saving plots
    """
    
    println("="^80)
    println("🔬 CTMRG Bond Dimension Convergence Benchmark (AKLT)")
    println("="^80)
    println("Parameters: a1=$(@sprintf("%.4f", a1)), a2=$(@sprintf("%.4f", a2))")
    println("Bond dimensions: $χ_values")
    println("Fixed iterations: $nsteps")
    println()
    
    # Create output directory
    if save_plots && !isdir(plot_dir)
        mkpath(plot_dir)
        println("📁 Created directory: $(plot_dir)")
    end
    
    # Storage for results
    free_energies = Float64[]
    computation_times = Float64[]
    
    println("🔄 Computing CTMRG for different bond dimensions...")
    
    for (i, χmax) in enumerate(χ_values)
        print("  χ = $χmax: ")
        
        try
            # Time the computation
            start_time = time()
            free_energy = -ctmrg_exact_FE_density(a1, a2; χmax=χmax, cutoff=cutoff, nsteps=nsteps)
            computation_time = time() - start_time
            
            push!(free_energies, free_energy)
            push!(computation_times, computation_time)
            
            println("f = $(@sprintf("%.8f", free_energy)), time = $(@sprintf("%.2f", computation_time))s")
            
        catch e
            println("❌ Error: $e")
            push!(free_energies, NaN)
            push!(computation_times, NaN)
        end
    end
    
    # Create plots
    println("\n📊 Creating bond dimension convergence plots...")
    
    # Plot 1: Free energy vs bond dimension
    p1 = plot(title="AKLT CTMRG: Free Energy vs Bond Dimension",
              xlabel="Bond Dimension χ", ylabel="Free Energy Density",
              xscale=:log10, legend=:bottomright, dpi=300)
    
    valid_mask = .!isnan.(free_energies)
    if any(valid_mask)
        plot!(p1, χ_values[valid_mask], free_energies[valid_mask],
              marker=:circle, markersize=6, linewidth=2,
              label="CTMRG ($(nsteps) steps)")
        
        # Add horizontal line for highest χ value as "converged" reference
        if sum(valid_mask) > 0
            converged_fe = free_energies[findlast(valid_mask)]
            hline!(p1, [converged_fe], linestyle=:dash, color=:red, alpha=0.7,
                   label="χ=$(χ_values[findlast(valid_mask)]) reference")
        end
    end
    
    # Plot 2: Computation time vs bond dimension
    p2 = plot(title="AKLT CTMRG: Computation Time vs Bond Dimension",
              xlabel="Bond Dimension χ", ylabel="Computation Time (s)",
              xscale=:log10, yscale=:log10, legend=:topleft, dpi=300)
    
    if any(valid_mask)
        plot!(p2, χ_values[valid_mask], computation_times[valid_mask],
              marker=:square, markersize=6, linewidth=2,
              label="CTMRG ($(nsteps) steps)", color=:blue)
    end
    
    # Plot 3: Convergence error vs bond dimension
    p3 = plot(title="AKLT CTMRG: Convergence Error vs Bond Dimension",
              xlabel="Bond Dimension χ", ylabel="Absolute Error",
              xscale=:log10, yscale=:log10, legend=:topright, dpi=300)
    
    if sum(valid_mask) > 1
        # Use highest χ as reference
        ref_fe = free_energies[findlast(valid_mask)]
        errors = abs.(free_energies .- ref_fe)
        
        # Remove the reference point itself
        error_mask = valid_mask .& (errors .> 1e-15)
        if any(error_mask)
            plot!(p3, χ_values[error_mask], errors[error_mask],
                  marker=:diamond, markersize=6, linewidth=2,
                  label="Error vs χ=$(χ_values[findlast(valid_mask)])", color=:red)
        end
    end
    
    # Display plots
    bond_dim_plot = plot(p1, p2, p3, layout=(3,1), size=(800, 900))
    display(bond_dim_plot)
    
    # Save if requested
    if save_plots
        savefig(bond_dim_plot, joinpath(plot_dir, "aklt_ctmrg_bond_dimension_convergence.png"))
        println("💾 Saved bond dimension plot to $(plot_dir)/")
    end
    
    return Dict(
        "chi_values" => χ_values,
        "free_energies" => free_energies,
        "computation_times" => computation_times,
        "a1" => a1,
        "a2" => a2,
        "nsteps" => nsteps
    )
end

function benchmark_ctmrg_iterations(;
                                   a1::Float64=sqrt(3/2), 
                                   a2::Float64=sqrt(6),
                                   χmax::Int=32,
                                   nsteps_values=[10, 25, 50, 100, 200, 400, 600, 800, 1000],
                                   cutoff::Float64=1e-14,
                                   save_plots::Bool=true,
                                   plot_dir::String="plots")
    """
    Benchmark CTMRG convergence vs iteration number for fixed bond dimension.
    
    Args:
        a1, a2: AKLT model parameters
        χmax: Fixed bond dimension
        nsteps_values: Numbers of iterations to test
        cutoff: Convergence cutoff
        save_plots: Whether to save plots
        plot_dir: Directory for saving plots
    """
    
    println("="^80)
    println("🔬 CTMRG Iteration Convergence Benchmark (AKLT)")
    println("="^80)
    println("Parameters: a1=$(@sprintf("%.4f", a1)), a2=$(@sprintf("%.4f", a2))")
    println("Fixed bond dimension: χ = $χmax")
    println("Iteration counts: $nsteps_values")
    println()
    
    # Create output directory
    if save_plots && !isdir(plot_dir)
        mkpath(plot_dir)
        println("📁 Created directory: $(plot_dir)")
    end
    
    # Storage for results
    free_energies = Float64[]
    computation_times = Float64[]
    
    println("🔄 Computing CTMRG for different iteration counts...")
    
    for (i, nsteps) in enumerate(nsteps_values)
        print("  $(nsteps) steps: ")
        
        try
            # Time the computation
            start_time = time()
            free_energy = -ctmrg_exact_FE_density(a1, a2; χmax=χmax, cutoff=cutoff, nsteps=nsteps)
            computation_time = time() - start_time
            
            push!(free_energies, free_energy)
            push!(computation_times, computation_time)
            
            println("f = $(@sprintf("%.8f", free_energy)), time = $(@sprintf("%.2f", computation_time))s")
            
        catch e
            println("❌ Error: $e")
            push!(free_energies, NaN)
            push!(computation_times, NaN)
        end
    end
    
    # Create plots
    println("\n📊 Creating iteration convergence plots...")
    
    # Plot 1: Free energy vs iterations
    p1 = plot(title="AKLT CTMRG: Free Energy vs Iterations",
              xlabel="Number of Iterations", ylabel="Free Energy Density",
              legend=:bottomright, dpi=300)
    
    valid_mask = .!isnan.(free_energies)
    if any(valid_mask)
        plot!(p1, nsteps_values[valid_mask], free_energies[valid_mask],
              marker=:circle, markersize=6, linewidth=2,
              label="CTMRG (χ=$χmax)")
        
        # Add horizontal line for highest iteration count as reference
        if sum(valid_mask) > 0
            converged_fe = free_energies[findlast(valid_mask)]
            hline!(p1, [converged_fe], linestyle=:dash, color=:red, alpha=0.7,
                   label="$(nsteps_values[findlast(valid_mask)]) steps reference")
        end
    end
    
    # Plot 2: Computation time vs iterations
    p2 = plot(title="AKLT CTMRG: Computation Time vs Iterations",
              xlabel="Number of Iterations", ylabel="Computation Time (s)",
              legend=:topleft, dpi=300)
    
    if any(valid_mask)
        plot!(p2, nsteps_values[valid_mask], computation_times[valid_mask],
              marker=:square, markersize=6, linewidth=2,
              label="CTMRG (χ=$χmax)", color=:blue)
    end
    
    # Plot 3: Convergence error vs iterations
    p3 = plot(title="AKLT CTMRG: Convergence Error vs Iterations",
              xlabel="Number of Iterations", ylabel="Absolute Error",
              yscale=:log10, legend=:topright, dpi=300)
    
    if sum(valid_mask) > 1
        # Use highest iteration count as reference
        ref_fe = free_energies[findlast(valid_mask)]
        errors = abs.(free_energies .- ref_fe)
        
        # Remove the reference point itself
        error_mask = valid_mask .& (errors .> 1e-15)
        if any(error_mask)
            plot!(p3, nsteps_values[error_mask], errors[error_mask],
                  marker=:diamond, markersize=6, linewidth=2,
                  label="Error vs $(nsteps_values[findlast(valid_mask)]) steps", color=:red)
        end
    end
    
    # Display plots
    iteration_plot = plot(p1, p2, p3, layout=(3,1), size=(800, 900))
    display(iteration_plot)
    
    # Save if requested
    if save_plots
        savefig(iteration_plot, joinpath(plot_dir, "aklt_ctmrg_iteration_convergence.png"))
        println("💾 Saved iteration plot to $(plot_dir)/")
    end
    
    return Dict(
        "nsteps_values" => nsteps_values,
        "free_energies" => free_energies,
        "computation_times" => computation_times,
        "a1" => a1,
        "a2" => a2,
        "chi_max" => χmax
    )
end

function full_ctmrg_benchmark(;
                             a1::Float64=sqrt(3/2), 
                             a2::Float64=sqrt(6),
                             save_plots::Bool=true,
                             plot_dir::String="plots")
    """
    Run complete CTMRG benchmark suite for AKLT model.
    """
    
    println("🚀 Running Complete CTMRG Benchmark Suite for AKLT Model")
    println("="^80)
    
    # Benchmark 1: Bond dimension convergence
    println("\n1️⃣ Bond Dimension Convergence Test")
    bond_results = benchmark_ctmrg_bond_dimension(
        a1=a1, a2=a2, save_plots=save_plots, plot_dir=plot_dir)
    
    # Benchmark 2: Iteration convergence  
    println("\n2️⃣ Iteration Convergence Test")
    iter_results = benchmark_ctmrg_iterations(
        a1=a1, a2=a2, save_plots=save_plots, plot_dir=plot_dir)
    
    # Summary plot combining key results
    println("\n📈 Creating summary comparison...")
    
    # Extract converged values for comparison
    bond_converged = bond_results["free_energies"][end]
    iter_converged = iter_results["free_energies"][end]
    
    println("\n📊 Summary Results:")
    println("="^50)
    println("Bond dimension convergence (χ=64, 400 steps): f = $(@sprintf("%.8f", bond_converged))")
    println("Iteration convergence (χ=32, 1000 steps): f = $(@sprintf("%.8f", iter_converged))")
    println("Difference: $(@sprintf("%.2e", abs(bond_converged - iter_converged)))")
    
    # Recommendations
    println("\n💡 Recommendations:")
    if !isnan(bond_converged) && !isnan(iter_converged)
        if abs(bond_converged - iter_converged) < 1e-6
            println("✅ Both approaches converged to similar values")
            println("   Recommended: χ=32, 400 steps for good accuracy/speed balance")
        else
            println("⚠️  Significant difference detected - may need higher χ or more steps")
        end
    end
    
    return Dict(
        "bond_dimension_results" => bond_results,
        "iteration_results" => iter_results,
        "parameters" => Dict("a1" => a1, "a2" => a2)
    )
end

full_ctmrg_benchmark(;a1=sqrt(3/2), 
                             a2=sqrt(6),
                             save_plots=true,
                             plot_dir="plots")