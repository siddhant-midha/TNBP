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
# (Assuming fig2.jl has similar structure to the other files)
# Add your specific parameters here

fig_dir = "fig2"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

# Data file path
data_file = joinpath(fig_dir, "data.jld2")

# Define parameters for computation and storage
βc = 0.362
βs = vcat([0.2, 0.3], collect(LinRange(0.34, 0.4, 10)), [0.4, 0.5, 0.6, 0.7, 0.8])
wts = [4,6,8,10]
wmax = maximum(wts)
alphas = LinRange(0.3,1.0,length(wts))
Ls = [20]

# Initialize storage for results
all_mean_contribs = Dict()

if load_only && isfile(data_file)
    # Load pre-computed data
    println("📂 Loading saved data from: $(data_file)")
    data = JLD2.load(data_file)
    # Extract data from the loaded file
    βs = data["βs"]
    wts = data["wts"] 
    Ls = data["Ls"]
    βc = data["βc"]
    alphas = data["alphas"]
    all_mean_contribs = data["all_mean_contribs"]
    println("✅ Data loaded successfully!")
else
    # Compute data from scratch
    println("🧮 Computing data...")
    
    for L in Ls
        N = 2 * L^2
        mean_contribs_by_L = Dict()
        
        cluster_data, cluster_filename = load_latest_single_site_cluster_file(; 
            size_filter="L$(L)", 
            weight_filter="w$(wmax)", 
            boundary_filter="periodic", 
            save_dir = "../../saved_clusters"
        )
        println("✅ Loaded cluster data: $(cluster_filename)")
        
        loop_objects = cluster_data.all_loops
        all_loops = [loop_object.edges for loop_object in loop_objects]

        # Precompute loop sets for each weight
        loops_by_weight = Dict{Int, Vector{Vector{Tuple{Int,Int}}}}()
        for wt in wts
            loops_by_weight[wt] = Vector{Vector{Tuple{Int,Int}}}()
        end
        for loop in all_loops
            loop_weight = length(loop)
            if haskey(loops_by_weight, loop_weight)
                push!(loops_by_weight[loop_weight], loop)
            end
        end

        mean_contribs_by_wt = Dict()
        
        for wt in wts
            mean_contribs = Float64[]
            for β in βs
                T = Ising2D.get_ising_tn(L, β)
                adj_mat, BPedges, links = BP.get_adj_mat(T)
                messages = BP.get_messages(T, BPedges, links; random_part=0.05)
                messages = BP.message_passing(T, messages, BPedges, adj_mat; α = .9, noise=0., max_iters=1000)
                Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
                T_normalized = BP.normalize_tensors(T, Z_l)

                vals = Float64[]
                for loop in loops_by_weight[wt]
                    loop_contri = abs(scalar(BP.loop_contribution(loop, messages, T_normalized, BPedges, links, adj_mat)))
                    push!(vals, loop_contri)
                end
                mean_val = isempty(vals) ? NaN : mean(vals)
                push!(mean_contribs, mean_val)
                
                # Show progress
                println("  Computed L=$(L), wt=$(wt), β=$(round(β, digits=3)): mean contrib = $(round(mean_val, digits=6))")
            end
            
            # Store results for this weight
            mean_contribs_by_wt[wt] = mean_contribs
        end
        
        # Store results for this L
        mean_contribs_by_L[L] = mean_contribs_by_wt
        all_mean_contribs[L] = mean_contribs_by_L[L]
    end
    
    # Save computed data
    println("💾 Saving computed data to: $(data_file)")
    JLD2.save(data_file, Dict(
        "βs" => βs,
        "wts" => wts,
        "Ls" => Ls,
        "βc" => βc,
        "alphas" => alphas,
        "all_mean_contribs" => all_mean_contribs
    ))
    println("✅ Data saved successfully!")
end

println("📊 Generating plots...")

# Generate plots using the loaded or computed data
for L in Ls
    p = plot(; 
        legend=:topright, 
        yscale=:log10, 
        fontfamily="Times New Roman",
        framestyle=:box, 
        grid=false, 
        thickness_scaling=1.5, 
        legendfont = font("Times New Roman", 8),
        guidefont = font("Times New Roman", 16),
        tickfont = font("Times New Roman", 8),
        size=(500, 500)
    )

    for (wi, wt) in enumerate(wts)
        # Get the mean contributions for this weight and L
        mean_contribs = all_mean_contribs[L][wt]
        
        plot!(p, βs, mean_contribs; 
            label="Loop weight $wt", 
            color=:steelblue, 
            alpha=alphas[wi], 
            linewidth=2, 
            fontfamily="Times New Roman"
        )
    end

    # Add vertical dashed line at critical beta (no legend entry)
    vline!(p, [βc]; linestyle=:dash, color=:black, label=false)

    savefig(p, joinpath(fig_dir, "mean_loop_contribution_vs_beta_L$(L).pdf"))
end

println("✅ Plots saved to: $(fig_dir)/ directory")
