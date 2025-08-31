include("../cluster_expansion_free_energy.jl")
include("../../functions/brute_force_new.jl")
using FileIO, Statistics



function get_cluster_contrb(L, max_weight, T_normalized, messages, edges, links, adj_mat)
    cluster_data, cluster_filename = load_latest_single_site_cluster_file(; size_filter="L$(L)", weight_filter="w$(max_weight)", boundary_filter="periodic", save_dir = "../../saved_clusters")

    cluster_contributions = compute_all_single_site_cluster_contributions(
        T_normalized, messages, edges, links, adj_mat, 
        cluster_data, L, max_weight)
    corrections = Dict{Int, ComplexF64}()
    
    for wt in [4,6,7,8,9,10]
        # println("🔍 Computing cluster expansion up to weight $wt")
        # Sum only contributions from clusters with weight <= max_weight
        correction = sum(contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= wt)
        corrections[wt] = correction
    end
    return corrections
end 

function loop_corrected_free_energy(L, T_normalized, messages, BPedges, BPlinks, BPadj_mat)
    N = 2 * L^2
    correction_l = Dict{Int, ComplexF64}()
    Z = 1.
    Z += ZCorrection4th(BP.loop_contribution,L,L,messages,T_normalized,BPedges,BPlinks,BPadj_mat)
    correction_l[4] = -log(Z)/N
    Z += ZCorrection6th(BP.loop_contribution,L,L,messages,T_normalized,BPedges,BPlinks,BPadj_mat)
    correction_l[6] = -log(Z)/N
    Z += ZCorrection7th(BP.loop_contribution,L,L,messages,T_normalized,BPedges,BPlinks,BPadj_mat)
    correction_l[7] = -log(Z)/N
    Z += ZCorrection8th(BP.loop_contribution,L,L,messages,T_normalized,BPedges,BPlinks,BPadj_mat)
    correction_l[8] = -log(Z)/N
    Z += ZCorrection9th(BP.loop_contribution,L,L,messages,T_normalized,BPedges,BPlinks,BPadj_mat)
    correction_l[9] = -log(Z)/N
    Z += ZCorrection10th(BP.loop_contribution,L,L,messages,T_normalized,BPedges,BPlinks,BPadj_mat)
    correction_l[10] = -log(Z)/N
    return correction_l
end
