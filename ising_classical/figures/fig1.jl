include("helpers.jl")
using Plots
using JLD2
using ArgParse
default(fontfamily="Times New Roman")

# Parse command-line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--load-only"
            help = "Load saved data instead of recomputing"
            action = :store_true
    end
    return parse_args(s)
end

args = parse_commandline()
load_only = args["load-only"]

# Parameters and setup
Ls = [20,30]  # Define Ls early so it can be used in file path
alphas = LinRange(0.3,1.0,length(Ls))
wts = [0,4,6,8,10]
max_weight = maximum(wts)
βs = [0.362]

fig_dir = "fig1"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

# Data file path - include L values in filename
data_file = joinpath(fig_dir, "data_L$(join(Ls, "_")).jld2")

# Define data structures to store results
all_results = Dict()

if load_only && isfile(data_file)
    # Load pre-computed data
    println("📂 Loading saved data from: $(data_file)")
    data = JLD2.load(data_file)
    # Extract your data from the loaded file
    Ls = data["Ls"]
    wts = data["wts"]
    βs = data["βs"]
    alphas = data["alphas"]
    all_cluster_errors = data["all_cluster_errors"]
    all_loop_errors = data["all_loop_errors"]
    println("✅ Data loaded successfully!")
else
    # Compute data from scratch
    println("🧮 Computing data...")
    
    # Store results for saving later
    all_cluster_errors = Dict()
    all_loop_errors = Dict()

    for β in βs
        cluster_errors_by_L = Dict()
        loop_errors_by_L = Dict()
        
        for (li, L) in enumerate(Ls)
            N = 2 * L^2
            T = Ising2D.get_ising_tn(L, β)
            adj_mat, BPedges, links = BP.get_adj_mat(T)
            messages = BP.get_messages(T, BPedges, links; random_part=0.05)
            messages = BP.message_passing(T, messages, BPedges, adj_mat;  α = .9, noise=0., max_iters=1000)
            println("✅ BP converged $(BP.check_self_consistency(T, messages, adj_mat))")
            Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
            T_normalized = BP.normalize_tensors(T, Z_l)
            bp_FE_density = -sum(log.(real.(Z_l))) / N
            exact_free_energy = Ising2D.free_energy(β)

            bp_error = abs(bp_FE_density - exact_free_energy)

            loop_errors = zeros(length(wts))
            cluster_errors = zeros(length(wts))

            loop_errors[1] = bp_error
            cluster_errors[1] = bp_error

            cluster_corrs = get_cluster_contrb(L, max_weight, T_normalized, messages, BPedges, links, adj_mat)
            loop_corrs = loop_corrected_free_energy(L, T_normalized, messages, BPedges, links, adj_mat)
            for (i,w) in enumerate(wts)
                if i > 1
                    loop_errors[i] = abs(bp_FE_density + loop_corrs[w] - exact_free_energy)
                    cluster_errors[i] = abs(bp_FE_density + cluster_corrs[w] - exact_free_energy)
                end 
            end
            
            # Store errors for this L value
            cluster_errors_by_L[L] = cluster_errors
            loop_errors_by_L[L] = loop_errors
        end
        
        # Store all errors for this β value
        all_cluster_errors[β] = cluster_errors_by_L
        all_loop_errors[β] = loop_errors_by_L
    end
    
    # Save computed data
    println("💾 Saving computed data to: $(data_file)")
    JLD2.save(data_file, Dict(
        "Ls" => Ls,
        "wts" => wts,
        "βs" => βs,
        "alphas" => alphas,
        "all_cluster_errors" => all_cluster_errors,
        "all_loop_errors" => all_loop_errors
    ))
    println("✅ Data saved successfully!")
end

println("📊 Generating plots...")

# Debug: Print data info
println("Weights: ", wts)
println("System sizes: ", Ls)
println("Beta values: ", βs)

# Generate plots using the loaded or computed data
for β in βs
    p = plot(;yscale=:log10, 
             legend=:bottomleft, 
             fontfamily="Times New Roman", 
             xticks=wts,  # Simplified tick specification
             framestyle=:box, 
             grid=false, 
             thickness_scaling=1.5, 
             legendfont = font("Times New Roman", 8),
             guidefont = font("Times New Roman", 16),
             tickfont = font("Times New Roman", 8),
             size=(500, 500)
            )
             
    for (li, L) in enumerate(Ls)
        # Get errors for this L and β
        cluster_errors = all_cluster_errors[β][L]
        loop_errors = all_loop_errors[β][L]
        
        # Debug: Print error info
        println("L=$L, cluster_errors: ", cluster_errors)
        println("L=$L, loop_errors: ", loop_errors)

        plot!(p, wts, cluster_errors;  
             label="Cluster L=$L", 
             marker=:circle, 
             color=:steelblue, 
             linewidth=2, 
             alpha=alphas[li])
             
        plot!(p, wts, loop_errors;  
             label="Loop L=$L", 
             marker=:square, 
             color=colorant"tomato", 
             linewidth=2, 
             alpha=alphas[li])
    end
    savefig(p, joinpath(fig_dir, "error_vs_weight_beta_$(round(β,digits=2)).pdf"))
end

println("✅ Plots saved to: $(fig_dir)/ directory")