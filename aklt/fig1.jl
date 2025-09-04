include("aklt.jl")
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("../functions/BP.jl")
using Serialization
using Plots

L = 20 
a1 = 0.5 
a2_range = vcat(collect(0.5:0.1:0.9), collect(0.9:0.01:1.1), collect(1.2:0.1:1.9), collect(1.9:0.002:2.1),collect(2.2:0.1:3.0)) 
                                                                     
max_weights = [4,6,8,10]
save_dir = "figs"
save_filename = "free_energy_vs_a2.pdf"
alphas = LinRange(0.3,1.0,length(max_weights[1:end]))
# Create save directory if it doesn't exist
if !isdir(save_dir)
    mkpath(save_dir)
end

save_path = joinpath(save_dir, save_filename)

println("📊 Creating free energy vs a2 plot using single-site clusters")
println("   Lattice size: $(L)×$(L)")
println("   Fixed a1: $a1")
println("   a2 range: $(a2_range[1]) to $(a2_range[end]) ($(length(a2_range)) points)")
println("   Weight truncations: $max_weights")

# Load single-site cluster data once
println("\n📂 Loading single-site cluster data once...")
global cluster_data = nothing
global cluster_filename = ""

try
    global cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L", weight_filter="w$(maximum(max_weights))")
    println("✅ Loaded single-site cluster data from: $(cluster_filename)")
catch e
    println("⚠️  Failed to load with specific weight, trying any weight...")
    try
        global cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L")
        println("✅ Loaded single-site cluster data from: $(cluster_filename)")
    catch e2
        error("❌ Could not load single-site cluster data for L=$L: $e2")
    end
end

# Verify cluster_data is not nothing
if cluster_data === nothing
    error("❌ Cluster data is still nothing after loading attempts")
end

# Initialize plot
p = plot(legend=:bottomleft, 
        fontfamily="Times New Roman",
        framestyle=:box, grid=false, thickness_scaling=1.5,
        legendfont = font("Times New Roman", 8),
        guidefont = font("Times New Roman", 16),
        tickfont = font("Times New Roman", 8),
        size=(500, 500))

# Compute exact CTMRG solution
println("\n🎯 Computing exact CTMRG solution...")
exact_f = [-ctmrg_exact_FE_density(a1, a2) for a2 in a2_range]
plot!(p, a2_range, exact_f, linewidth=3, color=:black, 
        label="CTMRG", linestyle=:solid)

# Initialize storage for all a2 points
bp_f = Float64[]
cluster_f_dict = Dict{Int, Vector{Float64}}()
for max_weight in max_weights
    cluster_f_dict[max_weight] = Float64[]
end

max_weight_overall = maximum(max_weights)

# Compute for each a2 value
println("\n📊 Computing for each a2 value...")
for (j, a2) in enumerate(a2_range)
    if j % 5 == 1 || j <= 3
        println("  Computing a2 = $a2 (point $j/$(length(a2_range)))...")
    end
    
    # Create AKLT tensor network
    T = aklt_norm_network(L; a1=a1, a2=a2)
    N = L^2
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    
    # Compute BP fixed point
    messages = BP.get_messages(T, edges, links; random_part=0.1)
    messages = BP.message_passing(T, messages, edges, adj_mat; α=1., max_iters=1000)
    
    # Get BP partition function (before normalization)
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Complex(Z_bp_full))
    
    # Normalize tensors
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)
    
    # BP free energy density
    f_bp = -real(log_Z_bp_full) / N
    push!(bp_f, f_bp)
    
    # Compute ALL cluster contributions once for this a2
    cluster_contributions = compute_all_single_site_cluster_contributions(
        T_normalized, messages, edges, links, adj_mat, 
        cluster_data, L, max_weight_overall)
    
    # Compute corrections for each weight
    for max_weight in max_weights
        correction = sum(contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= max_weight)
        f_corrected = f_bp + correction
        push!(cluster_f_dict[max_weight], f_corrected)
    end
end

# Plot BP baseline
plot!(p, a2_range, bp_f, linewidth=2, color=:steelblue, 
        label="BP Approximation", linestyle=:dash)

# Colors for different weights
colors = [:red, :blue, :green, :orange, :purple, :brown]

# Plot cluster expansion for each weight
for (i, max_weight) in enumerate(max_weights)
    color = colors[min(i, length(colors))]
    plot!(p, a2_range, cluster_f_dict[max_weight], linewidth=2, color=:steelblue,alpha=alphas[i],
            label="Cluster w≤$max_weight", linestyle=:solid)
end

# Save plot
if !isempty(save_path)
    savefig(p, save_path)
    println("✅ Main plot saved to: $save_path")
end

# Create error plot
println("\n📊 Creating error plot...")
p_error = plot(legend=:bottomright, 
        fontfamily="Times New Roman",
        framestyle=:box, grid=false, thickness_scaling=1.5,
        legendfont = font("Times New Roman", 8),
        guidefont = font("Times New Roman", 16),
        tickfont = font("Times New Roman", 8),
        size=(500, 500),
        yscale=:log10)

# Plot BP error
bp_error = abs.(bp_f .- exact_f)
# Filter out any NaN, Inf, or zero values that cause plotting issues
valid_indices = .!(isnan.(bp_error) .| isinf.(bp_error) .| (bp_error .<= 0))
if any(valid_indices)
    plot!(p_error, a2_range[valid_indices], bp_error[valid_indices], linewidth=2, color=:steelblue, 
            label="BP Error", linestyle=:dash)
end

# Plot cluster expansion errors for each weight
for (i, max_weight) in enumerate(max_weights)
    color = colors[min(i, length(colors))]
    cluster_error = abs.(cluster_f_dict[max_weight] .- exact_f)
    # Filter out any NaN, Inf, or zero values that cause plotting issues
    valid_indices = .!(isnan.(cluster_error) .| isinf.(cluster_error) .| (cluster_error .<= 0))
    if any(valid_indices)
        plot!(p_error, a2_range[valid_indices], cluster_error[valid_indices], linewidth=2, color=:steelblue,alpha=alphas[i],
                label="Cluster w≤$max_weight Error", linestyle=:solid)
    end
end

# Save error plot
if !isempty(save_path)
    error_save_path = replace(save_path, ".pdf" => "_errors.pdf")
    savefig(p_error, error_save_path)
    println("✅ Error plot saved to: $error_save_path")
end
