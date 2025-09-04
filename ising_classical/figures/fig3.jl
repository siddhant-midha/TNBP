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
# (Assuming fig3.jl has similar structure to the other files)
# Add your specific parameters here

fig_dir = "fig3"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

# Data file path
data_file = joinpath(fig_dir, "data.jld2")

# Define parameters
q = 3.9
βs = [0.2, 0.362, 0.5]
wts = [4,6,8,10]
wmax = maximum(wts)
Ls = [20]
color_grad = palette(:pastel, length(βs))

# Initialize storage for results
all_results = Dict()

if load_only && isfile(data_file)
    # Load pre-computed data
    println("📂 Loading saved data from: $(data_file)")
    data = JLD2.load(data_file)
    # Extract data
    βs = data["βs"]
    wts = data["wts"]
    Ls = data["Ls"]
    color_grad = data["color_grad"]
    all_results = data["all_results"]
    println("✅ Data loaded successfully!")
else
    # Compute data from scratch
    println("🧮 Computing data...")

    for L in Ls
        N = 2 * L^2
        results_by_L = Dict()
        
        println("Loading cluster data for L=$(L)...")
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
        
        results_by_beta = Dict()

        for (bi, β) in enumerate(βs)
            println("Computing for L=$(L), β=$(β)...")
            mean_contribs = Float64[]
            T = Ising2D.get_ising_tn(L, β)
            adj_mat, BPedges, links = BP.get_adj_mat(T)
            messages = BP.get_messages(T, BPedges, links; random_part=0.05)
            messages = BP.message_passing(T, messages, BPedges, adj_mat; α = .9, noise=0., max_iters=1000)
            Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
            T_normalized = BP.normalize_tensors(T, Z_l)

            # Store contributions by weight
            wt_contribs = Dict()
            
            for wt in wts
                vals = Float64[]
                for loop in loops_by_weight[wt]
                    loop_contri = abs(scalar(BP.loop_contribution(loop, messages, T_normalized, BPedges, links, adj_mat)))
                    push!(vals, loop_contri)
                end
                mean_val = isempty(vals) ? NaN : mean(vals)
                push!(mean_contribs, mean_val)
                
                # Store all individual values for this weight
                wt_contribs[wt] = vals
            end
            
            # Store means and all values for this beta
            results_by_beta[β] = Dict(
                "mean_contribs" => mean_contribs,
                "wt_contribs" => wt_contribs
            )
        end
        
        # Store results for this L
        results_by_L[L] = results_by_beta
        all_results[L] = results_by_L[L]
    end
    
    # Save computed data
    println("💾 Saving computed data to: $(data_file)")
    JLD2.save(data_file, Dict(
        "βs" => βs,
        "wts" => wts,
        "Ls" => Ls,
        "color_grad" => color_grad,
        "all_results" => all_results
    ))
    println("✅ Data saved successfully!")
end

println("📊 Generating plots...")

# Generate plots using the loaded or computed data
for L in Ls
    p = plot(; 
        yscale=:log10, 
        fontfamily="Times New Roman", 
        xticks=(wts, string.(wts)),
        framestyle=:box, grid=false, thickness_scaling=1.5,
        legendfont = font("Times New Roman", 8),
        guidefont = font("Times New Roman", 16),
        tickfont = font("Times New Roman", 8),
        size=(500, 500)
    )

    for (bi, β) in enumerate(βs)
        # Get the mean contributions for this beta and L
        mean_contribs = all_results[L][β]["mean_contribs"]
        
        plot!(p, wts, mean_contribs;  
              color=color_grad[bi], 
              marker=:circle, 
              linewidth=2, 
              fontfamily="Times New Roman", 
              label=L"\beta = %$(β)")
    end
    
    savefig(p, joinpath(fig_dir, "mean_loop_contribution_vs_weight_L$(L)_βmin$(round(minimum(βs),digits=2))_βmax$(round(maximum(βs),digits=2)).pdf"))
end

println("✅ Plots saved to: $(fig_dir)/ directory")
