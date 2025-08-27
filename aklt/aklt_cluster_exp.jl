using Plots
include("aklt.jl")
push!(LOAD_PATH, "../functions/")
using BP

## trivial deformation to get AKLT state 

a1 = sqrt(3/2)
a2 = sqrt(6)

L = 10
N = L^2

TN = aklt_norm_network(L; a1 = a1, a2 = a2)
adj_mat, edges, links = BP.get_adj_mat(TN)
messages = BP.get_messages(TN, edges, links; random_part=0.1)
messages = BP.message_passing(TN, messages, edges, adj_mat; α=0.8, max_iters=1000)
Z_list = BP.get_fixed_point_list(TN, messages, adj_mat)
TN = BP.normalize_tensors(TN, Z_list)
TN_normalized = BP.normalize_tensors(TN, Z_list)
Z_list_complex = Complex.(Z_list)
bp_FE_per_site = sum(log.(Z_list_complex)) / N
bp_PF_per_site = exp(bp_FE_per_site)
exact_FE_per_site = ctmrg_exact_FE_density(a1,a2)
bp_error = abs(bp_FE_per_site - exact_FE_per_site)

error_arr = [bp_error]
w_ar = [0,4,8,10]

for w in w_ar[2:end]
    cluster_data = load_latest_cluster_file(L, w)
    clustercorrx = cluster_contr_by_site(cluster_data, TN_normalized, messages, edges, links, adj_mat)
    cluster_FE_per_site_correction = clustercorrx / N 
    corrected_FE_per_site = bp_FE_per_site + cluster_FE_per_site_correction
    cluster_error = abs(corrected_FE_per_site - exact_FE_per_site)
    push!(error_arr, cluster_error)
end 

p = plot(w_ar, error_arr,
    marker=:circle,
    xlabel="wt",
    title="Error FE Density",
    legend=:topright,
    yscale=:log10,
    color=:red)
    
display(p)
readline()