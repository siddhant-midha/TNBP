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

# Parameters
L = 20  # Lattice size
wts = [0, 4, 6, 8, 10]  # Weights for cluster expansion (0 represents BP)
βs = collect(0.25:0.001:0.45)  # Focus on the phase transition region
q = 3.89
βc = 0.361 # 1/2 * log(q / (q-2))  # Critical beta
max_weight = maximum(wts)
alphas = LinRange(0.3,1.0,length(wts[2:end]))

# Replace the color setup with a more visually appealing scheme
# Colors for BP, and cluster weights 4,6,8,10
colors = [:steelblue, :darkgreen, :darkorange, :crimson, :indigo]  # Professional color scheme

# Create output directory
fig_dir = "fig4"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

# Data file path
data_file = joinpath(fig_dir, "data_L$(L).jld2")

# Initialize data arrays
N = 2 * L^2  # Number of sites
exact_FE_arr = zeros(length(βs))
bp_FE_arr = zeros(length(βs))
cluster_FE_arr = Dict(w => zeros(length(βs)) for w in wts[2:end])  # Skip 0 (BP)

if load_only && isfile(data_file)
    # Load pre-computed data
    println("📂 Loading saved data from: $(data_file)")
    data = JLD2.load(data_file)
    exact_FE_arr = data["exact_FE_arr"]
    bp_FE_arr = data["bp_FE_arr"]
    cluster_FE_arr = data["cluster_FE_arr"]
    println("✅ Data loaded successfully!")
else
    # Compute data from scratch
    println("🧮 Computing free energies across β range...")
    
    # Load cluster data for this system size
    cluster_data, cluster_filename = load_latest_single_site_cluster_file(
        size_filter="L$(L)", 
        weight_filter="w$(max_weight)", 
        boundary_filter="periodic", 
        save_dir="../../saved_clusters"
    )
    println("✅ Loaded cluster data: $(cluster_filename)")
    
    for (bi, β) in enumerate(βs)
        # Generate Ising tensor network
        T = Ising2D.get_ising_tn(L, β)
        adj_mat, BPedges, links = BP.get_adj_mat(T)
        
        # Perform BP
        messages = BP.get_messages(T, BPedges, links; random_part=0.05)
        messages = BP.message_passing(T, messages, BPedges, adj_mat; α=.9, noise=0., max_iters=1000)
        Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
        T_normalized = BP.normalize_tensors(T, Z_l)
        
        # Compute BP free energy
        bp_FE_density = -sum(log.(real.(Z_l))) / N
        bp_FE_arr[bi] = bp_FE_density
        
        # Compute exact free energy
        exact_FE = Ising2D.free_energy(β)
        exact_FE_arr[bi] = exact_FE
        
        # Compute cluster corrections for each weight
        cluster_corrs = get_cluster_contrb(L, max_weight, T_normalized, messages, BPedges, links, adj_mat)
        
        # Calculate cluster-corrected free energies
        for w in wts[2:end]  # Skip 0 (BP)
            cluster_FE_arr[w][bi] = bp_FE_density + cluster_corrs[w]
        end
        
        # Show progress
        if bi % 20 == 0 || bi == length(βs)
            println("  Progress: $(bi)/$(length(βs)) (β = $(β))")
        end
    end
    
    # Save computed data
    println("💾 Saving computed data to: $(data_file)")
    JLD2.save(data_file, Dict(
        "exact_FE_arr" => exact_FE_arr,
        "bp_FE_arr" => bp_FE_arr,
        "cluster_FE_arr" => cluster_FE_arr,
        "βs" => βs,
        "L" => L,
        "wts" => wts
    ))
    println("✅ Data saved successfully!")
end

println("📊 Generating plots...")

# Create the plots
# 1. Main plot: Free energies vs beta
p1 = plot(; 
    legend=:bottomleft, 
    fontfamily="Times New Roman",
    framestyle=:box, grid=false, thickness_scaling=1.5,
    legendfont = font("Times New Roman", 8),
    guidefont = font("Times New Roman", 16),
    tickfont = font("Times New Roman", 8),
    size=(500, 500) # Square plot with 1:1 aspect ratio
    )

# 2. Error plot: Absolute difference from exact solution
p2 = plot(; 
    legend=false, yscale=:log10,
    fontfamily="Times New Roman",
    framestyle=:box, grid=false, thickness_scaling=1.5,
    legendfont = font("Times New Roman", 8),
    guidefont = font("Times New Roman", 16),
    tickfont = font("Times New Roman", 8),
    size=(500, 500) # Square plot with 1:1 aspect ratio
    )

# Plot exact solution with a distinctive color and style
plot!(p1, βs, exact_FE_arr, label="Onsager", color=:black, linewidth=3, 
      fontfamily="Times New Roman", linestyle=:solid)

# Plot BP solution with a distinctive style
plot!(p1, βs, bp_FE_arr, label="BP vacuum", color=:steelblue, linewidth=2.5, linestyle=:dash,
      fontfamily="Times New Roman")

# Plot cluster solutions with distinct colors
for (i, w) in enumerate(wts[2:end])  # Skip 0 (BP)
    plot!(p1, βs, cluster_FE_arr[w], 
        label="Cluster weight $w", 
        color=:steelblue, 
        alpha=alphas[i],
        linewidth=2., 
        fontfamily="Times New Roman",
        linestyle=:solid
    )
end

# Add vertical line at critical beta
vline!(p1, [βc], linestyle=:dash, color=:grey, linewidth=1.5, label=false, alpha=0.75)

# Save main plot
savefig(p1, joinpath(fig_dir, "free_energy_comparison_L$(L).pdf"))

# Calculate and plot errors
bp_error = abs.(bp_FE_arr .- exact_FE_arr)
plot!(p2, βs, bp_error, label="BP vacuum", color=colors[1], linewidth=2.5, 
      fontfamily="Times New Roman", linestyle=:dash)

for (i, w) in enumerate(wts[2:end])  # Skip 0 (BP)
    cluster_error = abs.(cluster_FE_arr[w] .- exact_FE_arr)
    plot!(p2, βs, cluster_error, 
        label="Cluster weight $w", 
        color=:steelblue,
        alpha = alphas[i], 
        linewidth=2, 
        fontfamily="Times New Roman",
        linestyle=:solid
    )
end

# Add vertical line at critical beta in error plot
vline!(p2, [βc], linestyle=:dash, color=:grey, linewidth=1.5, label=false, alpha=0.75)

# Save error plot
savefig(p2, joinpath(fig_dir, "free_energy_error_L$(L).pdf"))

println("✅ Plots saved to: $(fig_dir)/ directory")
println("  - free_energy_comparison_L$(L).pdf")
println("  - free_energy_error_L$(L).pdf")
