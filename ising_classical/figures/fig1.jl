include("helpers.jl")


Ls = [8, 10]
alphas = LinRange(0.3,1.0,length(Ls))
wts = [0,4,6,7,8,9,10]
max_weight = maximum(wts)


βs = [0.1,0.2,0.3,0.364,0.4,0.5,0.6,0.7,0.8]
fig_dir = "fig1"

if !isdir(fig_dir)
    mkpath(fig_dir)
end

for β in βs
    p = plot(; xlabel="Max cluster/loop weight", ylabel="Error", yscale=:log10, legend=:topright, title="BP + Corrections Error vs Weight, β=$(β)")
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

        plot!(p, wts, cluster_errors; label="Cluster L=$L", marker=:circle, color=:blue, linewidth=2, alpha=alphas[li])
        plot!(p, wts, loop_errors; label="Loop L=$L", marker=:square, color=:red, linewidth=2, alpha=alphas[li])
    end
    savefig(p, joinpath(fig_dir, "error_vs_weight_beta_$(β).pdf"))
end