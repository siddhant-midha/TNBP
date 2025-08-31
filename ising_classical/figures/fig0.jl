include("helpers.jl")

Ls = [10]
wts = [0,4,6,8,10]
color_wts = palette(:tab10, length(wts))  # Use a distinct color for each cluster order
alphas = LinRange(0.2,1.0,length(wts))

max_weight = maximum(wts)

βs = collect(0.2:0.001:0.6)
fig_dir = "fig0"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

for (li, L) in enumerate(Ls)
    N = 2 * L^2
    # Store free energies for each cluster order and exact
    cluster_FE_arr = zeros(length(wts), length(βs))
    exact_FE_arr = zeros(length(βs))
    # Compute free energies
    for (bi, β) in enumerate(βs)
        T = Ising2D.get_ising_tn(L, β)
        adj_mat, BPedges, links = BP.get_adj_mat(T)
        messages = BP.get_messages(T, BPedges, links; random_part=0.05)
        messages = BP.message_passing(T, messages, BPedges, adj_mat;  α = .9, noise=0., max_iters=1000)
        Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
        T_normalized = BP.normalize_tensors(T, Z_l)
        bp_FE_density = -sum(log.(real.(Z_l))) / N
        exact_FE = Ising2D.free_energy(β) #/ β
        exact_FE_arr[bi] = exact_FE

        cluster_corrs = get_cluster_contrb(L, max_weight, T_normalized, messages, BPedges, links, adj_mat)
        cluster_FE_arr[1, bi] = bp_FE_density  #/ β
        for (i, w) in enumerate(wts)
            if i > 1
                cluster_FE_arr[i, bi] = (bp_FE_density + cluster_corrs[w])  #/ β
            end
        end
    end

    # Plot free energies vs beta
    p1 = plot(; xlabel="β", ylabel="Free Energy", title="Free Energy vs β, L=$L", legend=:topright)
    plot!(p1, βs, exact_FE_arr; label="Exact", color=:black, linewidth=2)
    for (i, w) in enumerate(wts)
        plot!(p1, βs, cluster_FE_arr[i, :]; label="Cluster order $w", color=color_wts[i], linewidth=2)
    end
    savefig(p1, joinpath(fig_dir, "free_energy_vs_beta_L$(L).png"))

    # Compute second derivative w.r.t. beta (numerical)
    d2_exact = [NaN; diff(diff(exact_FE_arr)) ./ (diff(βs)[1])^2]

    # Compute second derivatives for cluster orders, avoiding reshape error
    d2_cluster = zeros(length(wts), length(βs))
    for i in 1:length(wts)
        # Second derivative for cluster order i
        d2_vals = [NaN; diff(diff(cluster_FE_arr[i, :])) ./ (diff(βs)[1])^2]
        # Pad to length(βs) with NaN if needed
        if length(d2_vals) < length(βs)
            d2_vals = vcat(d2_vals, fill(NaN, length(βs) - length(d2_vals)))
        end
        d2_cluster[i, :] .= d2_vals
    end

    # Plot second derivative vs beta
    p2 = plot(; xlabel="β", ylabel="d²F/dβ²", title="Second Derivative of Free Energy vs β, L=$L", legend=:topright)
    plot!(p2, βs, d2_exact; label="Exact", color=:black, linewidth=2)
    for (i, w) in enumerate(wts)
        plot!(p2, βs, d2_cluster[i, :]; label="Cluster order $w", color=color_wts[i], linewidth=2)
    end
    savefig(p2, joinpath(fig_dir, "second_derivative_free_energy_vs_beta_L$(L).pdf"))
end