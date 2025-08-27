using Plots
include("aklt.jl")
push!(LOAD_PATH, "../functions/")
using BP


L = 10
N = L^2

w = 10
cluster_data = load_latest_cluster_file(L, w)
all_loops = cluster_data.all_loops;  
all_loops_ = [loop_object.edges for loop_object in all_loops]; # this contains list of loops as list of edges 
clusters_by_site = cluster_data.clusters_by_site;

a1 = 0.5
a2_arr = 0.8:0.01:1.2

loop_arr = []
error_FE_arr = []

@showprogress for a2 in a2_arr 
    TN = aklt_norm_network(L; a1 = a1, a2 = a2)
    adj_mat, edges, links = BP.get_adj_mat(TN)
    messages = BP.get_messages(TN, edges, links; random_part=0.1)
    messages = BP.message_passing(TN, messages, edges, adj_mat; α=0.75, max_iters=1000)
    Z_list = BP.get_fixed_point_list(TN, messages, adj_mat)
    TN = BP.normalize_tensors(TN, Z_list)
    TN_normalized = BP.normalize_tensors(TN, Z_list)
    Z_list_complex = Complex.(Z_list)
    bp_FE_per_site = sum(log.(Z_list_complex)) / N
    bp_PF_per_site = exp(bp_FE_per_site)
    exact_FE_per_site = ctmrg_exact_FE_density(a1,a2)
    loopsum = abs(sum([scalar(BP.loop_contribution(loop, messages, TN_normalized, edges, links, adj_mat)) for loop in all_loops_]))
    push!(loop_arr, loopsum)
    push!(error_FE_arr,abs(exact_FE_per_site - bp_FE_per_site) )
end 

# Create two plots sharing x-axis
p1 = plot(a2_arr, loop_arr,
    marker=:square,
    label="Loop Sum",
    xlabel="a2",
    ylabel="Loop Sum",
    legend=:topright,
    yscale=:ln,
    color=:blue)

p2 = plot(a2_arr, error_FE_arr,
    marker=:circle,
    label="Error FE",
    xlabel="a2",
    ylabel="Error FE",
    legend=:topright,
    yscale=:ln,
    color=:red)

# Combine plots vertically with linked x-axis
p = plot(p1, p2, layout=(2,1), link=:x, size=(700, 500))

display(p)