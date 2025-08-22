include("randompeps.jl")

# Ursell function for connected clusters (from original implementation)
function ursell_function(cluster::Cluster)
    """
    Compute the Ursell function φ(W) for a connected cluster W.
    For connected clusters, this is (-1)^(|W|-1) * (|W|-1)!
    where |W| is the total number of loops in the cluster (with multiplicities).
    """
    total_loops = cluster.total_loops
    if total_loops == 1
        return 1.0
    else
        # φ(W) = (-1)^(|W|-1) * (|W|-1)!
        sign = (total_loops - 1) % 2 == 0 ? 1.0 : -1.0
        factorial_part = factorial(big(total_loops - 1))
        return sign * Float64(factorial_part)
    end
end

function cluster_contr(T_normalized, messages, edges, links, adj_mat, cluster, all_loops)
    # Compute Ursell function φ(W)
    phi_W = ursell_function(cluster)
    
    # if abs(phi_W) < 1e-15
    #     return 0.0  # Skip negligible contributions
    # end
    
    # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
    Z_W = 1.0 + 0im  # Complex number for cluster contribution
    # invalidloop = false 
    for loop_id in cluster.loop_ids
        multiplicity = cluster.multiplicities[loop_id]
        loop = all_loops[loop_id]
        
        # Convert loop edges format for BP.jl (ensure v1 < v2 ordering)
        loop_edges_bp = Tuple{Int,Int}[]
        for edge in loop.edges
            v1, v2 = edge
            push!(loop_edges_bp, (min(v1, v2), max(v1, v2)))
        end
        
        # Compute loop contribution Z_l using BP
        try
            Z_l_tensor = BP.loop_contribution(loop_edges_bp, messages, T_normalized, edges, links, adj_mat)
            # Extract scalar value from ITensor
            Z_l = scalar(Z_l_tensor)
            Z_W *= Z_l^multiplicity
        catch e
            # println("⚠️  Error computing loop contribution for loop $loop_id: $e")
            # println("invalid loops present...", !(all([e in edges for e in loop_edges_bp])))
            return 0.0  # Return 0 instead of breaking
        end
    end
    
    # Add contribution to free energy density
    contribution = -phi_W * Z_W 
    return real(contribution)  # Return real part only
end 


# Fixed parameters
ti = true   
orthog = false   
normalise = true 
annealing = .75 
maxiter = 1000


# w_list = [4,6,8,10]
# # Load loop data for all weights
# cluster_data_diff_weights = Dict()
# for w in w_list
#     data = load_latest_cluster_file(N,w)
#     cluster_data_diff_weights[w] = data
# end


## Let me do this for a fixed weight w, fixed N, fixed \eta, just one sample of rPEPS
N = 10
T = N
η = 0.1
nsamples = 1
w = 4

cluster_data = load_latest_cluster_file(N,w)
clusters_by_site = cluster_data.clusters_by_site
all_loops = cluster_data.all_loops

tensors, peps = peps_controllable(N, T; η=η, ti=ti, orthog=orthog)
exact_fe = log(real(contract_peps_no_phys(peps; cutoff=1E-12, maxdim=2^N)))

adj_mat, edges, links = BP.get_adj_mat(tensors)
messages = BP.get_messages(tensors, edges, links; random_part=0.0)
messages = BP.message_passing(tensors, messages, edges, adj_mat; α = annealing, max_iters=maxiter)
println("self consistency chec...", BP.check_self_consistency(tensors, messages, adj_mat))
Z_l = BP.get_fixed_point_list(tensors, messages, adj_mat)
log_Z_bp_full = sum(log.(real.(Z_l)))
T_normalized = BP.normalize_tensors(tensors, Z_l)


global clustercorrx = 0
for site in 1:(N*T)
    clusters = clusters_by_site[site]
    for cluster in clusters
        contribution = cluster_contr(T_normalized, messages, edges, links, adj_mat, cluster, all_loops)
        # Only add contribution if computation was successful (non-zero)
        if !isnan(contribution) && isfinite(contribution)
            global clustercorrx += contribution
        end
    end
end 

println("BP error: ", abs(log_Z_bp_full - exact_fe))
println("BP + cluster correction error: ", abs(log_Z_bp_full + clustercorrx - exact_fe))