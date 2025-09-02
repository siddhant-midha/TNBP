include("helpers.jl")
using Plots
default(fontfamily="Times New Roman")

Ls = [10]
wts = [0,4,6,8,10]
color_wts = palette(:tab10, length(wts))  # Use a distinct color for each cluster order
alphas = LinRange(0.2,1.0,length(wts))

max_weight = maximum(wts)

βs = collect(0.1:0.001:0.6)
fig_dir = "fig0"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

for (li, L) in enumerate(Ls)
    N = 2 * L^2
    cluster_FE_arr = zeros(length(wts), length(βs))
    exact_FE_arr = zeros(length(βs))
    bp_FE_arr = zeros(length(βs))
    abs_error_arr = zeros(length(βs))
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
    end

    # Find the beta value where the error is maximum
    max_error_idx = argmax(abs_error_arr)
    max_error_beta = βs[max_error_idx]
    max_error_value = abs_error_arr[max_error_idx]
    
    println("Maximum error occurs at β = $(max_error_beta)")
    println("Maximum error value: $(max_error_value)")
    
    # Plot free energies and absolute error vs beta (dual y-axis)
    # Free energies on left y-axis, error on right y-axis
    p1 = plot(βs, exact_FE_arr, 
             label="Onsager", color=:red, linewidth=2, 
             fontfamily="Times New Roman", legend=:bottomleft)
    plot!(p1, βs, bp_FE_arr, 
         label="BP vacuum", color=:blue, linewidth=2, 
         fontfamily="Times New Roman", ylabel="Free Energy")
    
    # Create a second y-axis for error (right side)
    p1_right = twinx(p1) # Create the right y-axis
    plot!(p1_right, βs, abs_error_arr, color=:grey, linewidth=2, 
          fontfamily="Times New Roman", linestyle = :dash, legend=false)
    
    # Add vertical line at maximum error
    # vline!(p1, [max_error_beta], color=:black, linestyle=:dot, label="Max error β=$(round(max_error_beta, digits=4))")
    # vline!(p1_right, [max_error_beta], color=:black, linestyle=:dot, label="")
    
    savefig(p1, joinpath(fig_dir, "free_energy_and_error_vs_beta_L$(L).pdf"))

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

    # Plot second derivative vs beta
    p2 = plot(; legend=:topright, fontfamily="Times New Roman")
    plot!(p2, βs, d2_exact; label="Onsager", color=:black, linewidth=2, fontfamily="Times New Roman")
    plot!(p2, βs, d2_bp; label="BP vacuum", color=:blue, linewidth=2, fontfamily="Times New Roman")
    
    # Add vertical line at maximum error to second plot
    vline!(p2, [max_error_beta], color=:black, linestyle=:dot, label="Max error β=$(round(max_error_beta, digits=4))")
    
    savefig(p2, joinpath(fig_dir, "second_derivative_free_energy_vs_beta_L$(L).pdf"))
end