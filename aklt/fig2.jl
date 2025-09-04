include("aklt.jl")
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("../functions/BP.jl")
using Serialization
using Plots

L = 20 
# List of (a1, a2) pairs to analyze
# parameter_pairs = [(0.5, 0.5), (0.5, 1.5), (0.5, 3.0)] ## in phase 
parameter_pairs = [(0.5, 1.), (0.5, 2.0)] ## critical points 

max_weights = [0,4,6,8,10]
save_dir = "figs"
save_filename = "error_vs_weight_multiple_params.pdf"

# Create save directory if it doesn't exist
if !isdir(save_dir)
    mkpath(save_dir)
end

save_path = joinpath(save_dir, save_filename)
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

# Define pastel colors for different parameter pairs
pastel_colors = [:lightblue, :lightpink, :lightgreen, :lightyellow, :lightcoral, :lightgray]

# Create plot
println("📊 Creating plot...")
p = plot(legend=(0.6, 0.5), 
        fontfamily="Times New Roman",
        framestyle=:box, grid=false, thickness_scaling=1.5,
        legendfont = font("Times New Roman", 8),
        guidefont = font("Times New Roman", 16),
        tickfont = font("Times New Roman", 8),
        size=(500, 500),
        yscale=:log10,
        xticks=(max_weights, string.(max_weights)))

# Collect all errors to determine y-axis range
all_errors_for_range = Float64[]

# Process each parameter pair
for (idx, (a1, a2)) in enumerate(parameter_pairs)
    println("📊 Processing (a1=$a1, a2=$a2)...")
    
    # Create AKLT tensor network
    T = aklt_norm_network(L; a1=a1, a2=a2)
    N = L^2

    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)

    # Compute BP fixed point
    messages = BP.get_messages(T, edges, links; random_part=0.1)
    messages = BP.message_passing(T, messages, edges, adj_mat; α=0.8, max_iters=1000)

    # Get BP partition function and free energy density
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Complex(Z_bp_full))
    f_bp = -real(log_Z_bp_full) / N

    # Normalize tensors
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)

    # Compute exact CTMRG solution
    exact_fe = -ctmrg_exact_FE_density(a1, a2)

    # Compute cluster contributions
    cluster_contributions = compute_all_single_site_cluster_contributions(
        T_normalized, messages, edges, links, adj_mat, 
        cluster_data, L, maximum(max_weights))

    # Compute errors for each weight
    weights = [0; max_weights]  # Include BP (weight 0)
    errors = zeros(length(weights))

    # BP error (weight 0)
    errors[1] = abs(f_bp - exact_fe)

    # Cluster expansion errors
    for (i, w) in enumerate(max_weights)
        if w == 0
            # For weight 0, no cluster correction (same as BP)
            errors[i+1] = errors[1]
        else
            # Filter contributions by weight and sum, with default 0 if empty
            matching_contribs = [contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= w]
            correction = isempty(matching_contribs) ? 0.0 : sum(matching_contribs)
            f_corrected = f_bp + correction
            errors[i+1] = abs(f_corrected - exact_fe)
        end
    end

    # Collect errors for range calculation
    append!(all_errors_for_range, errors[errors .> 0])

    # Plot errors for this parameter pair
    color = pastel_colors[min(idx, length(pastel_colors))]
    plot!(p, weights, errors, linewidth=2, color=color, 
            marker=:circle, markersize=4, 
            label=L"(a_1=%$a1, a_2=%$a2)")
end

# Set y-axis limits to powers of 10 based on all data
min_error = minimum(all_errors_for_range)
max_error = maximum(all_errors_for_range)
y_min_power = floor(Int, log10(min_error))
y_max_power = ceil(Int, log10(max_error))
y_ticks = 10.0 .^ (y_min_power:y_max_power)

plot!(p, yticks=y_ticks, ylims=(10.0^y_min_power, 10.0^y_max_power))

# Save plot
savefig(p, save_path)
println("✅ Plot saved to: $save_path")
