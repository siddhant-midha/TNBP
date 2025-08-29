include("../../functions/ClusterEnumeration.jl")
include("../../functions/Ising2D.jl")
include("../../functions/BP.jl")
include("../../functions/brute_force_new.jl")

using Colors, Plots, Graphs, ITensors, Statistics

βs = vcat([0.2, 0.3], collect(LinRange(0.34, 0.4, 10)), [0.4, 0.44, 0.5, 0.6, 0.7, 0.8])

wts = [4,6,8,10]
wmax = maximum(wts)
color_grad = palette(:viridis, length(wts))

files = Dict{Int,String}()
files[6] = "periodic_clusters_L6_w$(wmax).jld2" 
files[8] = "periodic_clusters_L8_w$(wmax).jld2" 
files[10] = "periodic_clusters_L10_w$(wmax).jld2" 


fig_dir = "figs"
if !isdir(fig_dir)
    mkpath(fig_dir)
end

Ls = [6,8,10]
for L in Ls
    N = 2 * L^2
    p = plot(; xlabel="β", ylabel="Mean loop contribution", title="Mean Loop Contribution vs β for different loop orders, L=$L", legend=:topright, yscale=:log10)
    loaded_data = open(files[L], "r") do io
        deserialize(io)
    end
    cluster_data = loaded_data["data"]
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

    for (wi, wt) in enumerate(wts)
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
        end
        plot!(p, βs, mean_contribs; label="Loop order $wt", color=color_grad[wi], marker=:circle, linewidth=2)
    end
    savefig(p, joinpath(fig_dir, "mean_loop_contribution_vs_beta_L$(L).png"))
end
