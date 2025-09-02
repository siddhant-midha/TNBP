include("helpers.jl")
using Plots
default(fontfamily="Times New Roman")

# Parameters
L = 10  # Lattice size
wts = [0, 4, 6, 8, 10]  # Weights for cluster expansion (0 represents BP)
βs = collect(0.25:0.001:0.45)  # Focus on the phase transition region
q = 3.89
βc = 1/2 * log(q / (q-2))  # Critical beta
max_weight = maximum(wts)

# Replace the color setup with a more visually appealing scheme
# Colors for BP, and cluster weights 4,6,8,10
colors = [:steelblue, :darkgreen, :darkorange, :crimson, :indigo]  # Professional color scheme
# Or alternatively use ColorSchemes package for better gradients
# using ColorSchemes
# colors = [colorscheme_palette(:cool, 5)...]

# Create output directory
fig_dir = "fig4"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

# Create the plots
# 1. Main plot: Free energies vs beta
p1 = plot(; 
    legend=:bottomleft, 
    fontfamily="Times New Roman"
)

# 2. Error plot: Absolute difference from exact solution
p2 = plot(; 
    legend=:bottomleft, yscale=:log10,
    fontfamily="Times New Roman"
)

# Compute free energies and errors
N = 2 * L^2  # Number of sites
exact_FE_arr = zeros(length(βs))
bp_FE_arr = zeros(length(βs))
cluster_FE_arr = Dict(w => zeros(length(βs)) for w in wts[2:end])  # Skip 0 (BP)

# Load cluster data for this system size
cluster_data, cluster_filename = load_latest_single_site_cluster_file(
    size_filter="L$(L)", 
    weight_filter="w$(max_weight)", 
    boundary_filter="periodic", 
    save_dir="../../saved_clusters"
)
println("✅ Loaded cluster data: $(cluster_filename)")

println("🧮 Computing free energies across β range...")
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

# Plot exact solution with a distinctive color and style
plot!(p1, βs, exact_FE_arr, label="Exact", color=:black, linewidth=3, 
      fontfamily="Times New Roman", linestyle=:solid)

# Plot BP solution with a distinctive style
plot!(p1, βs, bp_FE_arr, label="BP", color=colors[1], linewidth=2.5, 
      fontfamily="Times New Roman", linestyle=:dash)

# Plot cluster solutions with distinct colors
for (i, w) in enumerate(wts[2:end])  # Skip 0 (BP)
    plot!(p1, βs, cluster_FE_arr[w], 
        label="Cluster w≤$w", 
        color=colors[i+1], 
        linewidth=2, 
        fontfamily="Times New Roman",
        linestyle=:solid
    )
end

# Add vertical line at critical beta
vline!(p1, [βc], linestyle=:dash, color=:black, linewidth=1.5, label=false, alpha=0.75)

# Save main plot
savefig(p1, joinpath(fig_dir, "free_energy_comparison_L$(L).pdf"))

# Calculate and plot errors
bp_error = abs.(bp_FE_arr .- exact_FE_arr)
plot!(p2, βs, bp_error, label="BP", color=colors[1], linewidth=2.5, 
      fontfamily="Times New Roman", linestyle=:dash)

for (i, w) in enumerate(wts[2:end])  # Skip 0 (BP)
    cluster_error = abs.(cluster_FE_arr[w] .- exact_FE_arr)
    plot!(p2, βs, cluster_error, 
        label="Cluster w≤$w", 
        color=colors[i+1], 
        linewidth=2, 
        fontfamily="Times New Roman",
        linestyle=:solid
    )
end

# # Add vertical line at critical beta in error plot
vline!(p2, [βc], linestyle=:dash, color=:black, linewidth=1.5, label=false, alpha=0.75)

# Save error plot
savefig(p2, joinpath(fig_dir, "free_energy_error_L$(L).pdf"))

println("✅ Plots saved to: $(fig_dir)/ directory")
println("  - free_energy_comparison_L$(L).pdf")
println("  - free_energy_error_L$(L).pdf")
