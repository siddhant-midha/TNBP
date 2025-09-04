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

Ls = [20]
wts = [0,4,6,8,10]
color_wts = palette(:tab10, length(wts))  # Use a distinct color for each cluster order
alphas = LinRange(0.2,1.0,length(wts))

max_weight = maximum(wts)

βs = collect(0.1:0.01:0.6)
fig_dir = "fig0"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

# Data file path - include L values in filename
data_file = joinpath(fig_dir, "data_L$(join(Ls, "_")).jld2")

# Initialize data structures to store results
all_results = Dict()

for (li, L) in enumerate(Ls)
    N = 2 * L^2
    cluster_FE_arr = zeros(length(wts), length(βs))
    exact_FE_arr = zeros(length(βs))
    bp_FE_arr = zeros(length(βs))
    abs_error_arr = zeros(length(βs))
    
    if load_only && isfile(data_file)
        # Load pre-computed data
        println("📂 Loading saved data from: $(data_file)")
        data = JLD2.load(data_file)
        exact_FE_arr = data["exact_FE_arr_L$(L)"]
        bp_FE_arr = data["bp_FE_arr_L$(L)"]
        abs_error_arr = data["abs_error_arr_L$(L)"]
        println("✅ Data loaded successfully!")
    else
        # Compute data from scratch
        println("🧮 Computing free energies across β range for L=$(L)...")
        
        for (bi, β) in enumerate(βs)
            T = Ising2D.get_ising_tn(L, β)
            adj_mat, BPedges, links = BP.get_adj_mat(T)
            messages = BP.get_messages(T, BPedges, links; random_part=0.05)
            messages = BP.message_passing(T, messages, BPedges, adj_mat;  α = .9, noise=0., max_iters=1000)
            Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
            T_normalized = BP.normalize_tensors(T, Z_l)
            bp_FE_density = -sum(log.(real.(Z_l))) / N
            exact_FE = Ising2D.free_energy(β)
            bp_FE_arr[bi] = bp_FE_density
            exact_FE_arr[bi] = exact_FE
            abs_error_arr[bi] = abs(bp_FE_density - exact_FE)
            
            # Show progress
            if bi % 50 == 0 || bi == length(βs)
                println("  Progress: $(bi)/$(length(βs)) (β = $(β))")
            end
        end
        
        # Store results for this L value
        all_results["exact_FE_arr_L$(L)"] = exact_FE_arr
        all_results["bp_FE_arr_L$(L)"] = bp_FE_arr
        all_results["abs_error_arr_L$(L)"] = abs_error_arr
    end
    
    # Compute second derivative w.r.t. beta (numerical)
    d2_exact = [NaN; diff(diff(exact_FE_arr)) ./ (diff(βs)[1])^2]
    d2_bp = [NaN; diff(diff(bp_FE_arr)) ./ (diff(βs)[1])^2]
    # Pad to length(βs) with NaN if needed
    if length(d2_exact) < length(βs)
        d2_exact = vcat(d2_exact, fill(NaN, length(βs) - length(d2_exact)))
    end
    if length(d2_bp) < length(βs)
        d2_bp = vcat(d2_bp, fill(NaN, length(βs) - length(d2_bp)))
    end

    # --- find beta where |d2_bp| is maximized and print it ---
    tmp_abs = [isnan(x) ? -Inf : abs(x) for x in d2_bp]
    idx_max_abs = argmax(tmp_abs)
    beta_d2bp_max = βs[idx_max_abs]
    max_d2bp_value = d2_bp[idx_max_abs]
    println("→ Beta of max |d2_bp|: $(beta_d2bp_max)  (d2_bp = $(max_d2bp_value))")
    
    all_results["d2_exact_L$(L)"] = d2_exact
    all_results["d2_bp_L$(L)"] = d2_bp
    all_results["beta_d2bp_max_L$(L)"] = beta_d2bp_max
end

# Save computed data if not in load-only mode
if !load_only || !isfile(data_file)
    println("💾 Saving computed data to: $(data_file)")
    JLD2.save(data_file, all_results)
    println("✅ Data saved successfully!")
end

println("📊 Generating plots...")

for (li, L) in enumerate(Ls)
    # Get data for this L
    exact_FE_arr = load_only ? JLD2.load(data_file)["exact_FE_arr_L$(L)"] : all_results["exact_FE_arr_L$(L)"]
    bp_FE_arr = load_only ? JLD2.load(data_file)["bp_FE_arr_L$(L)"] : all_results["bp_FE_arr_L$(L)"]
    d2_exact = load_only ? JLD2.load(data_file)["d2_exact_L$(L)"] : all_results["d2_exact_L$(L)"]
    d2_bp = load_only ? JLD2.load(data_file)["d2_bp_L$(L)"] : all_results["d2_bp_L$(L)"]
    beta_d2bp_max = load_only ? JLD2.load(data_file)["beta_d2bp_max_L$(L)"] : all_results["beta_d2bp_max_L$(L)"]
    
    # Plot free energies and absolute error vs beta
    p1 = plot(βs, exact_FE_arr, 
             label="Onsager", color=:black, linewidth=2, 
             fontfamily="Times New Roman", legend=:bottomleft,
             framestyle=:box, grid=false, thickness_scaling=1.5,
             legendfont = font("Times New Roman", 8),
             guidefont = font("Times New Roman", 16),
             tickfont = font("Times New Roman", 8),
             size=(500, 500))
    plot!(p1, βs, bp_FE_arr, 
         label="BP vacuum", color=:steelblue, linewidth=2, 
         fontfamily="Times New Roman")
    
    savefig(p1, joinpath(fig_dir, "free_energy_and_error_vs_beta_L$(L).pdf"))

    # Plot second derivative vs beta
    p2 = plot(; legend=:topright,
             framestyle=:box, grid=false, thickness_scaling=1.5,
             legendfont = font("Times New Roman", 7),
             guidefont = font("Times New Roman", 7),
             tickfont = font("Times New Roman", 6),
             size=(500, 500))
    plot!(p2, βs, d2_exact; label="Onsager", color=:black, linewidth=2, fontfamily="Times New Roman")
    plot!(p2, βs, d2_bp; label="BP vacuum", color=:steelblue, linewidth=2, fontfamily="Times New Roman")
    
    # Add vertical dashed line at beta where |d2_bp| is maximal (no legend entry)
    vline!(p2, [beta_d2bp_max]; color=:grey, linestyle=:dash, label=false)

    savefig(p2, joinpath(fig_dir, "second_derivative_free_energy_vs_beta_L$(L).pdf"))
end

println("✅ Plots saved to: $(fig_dir)/ directory")