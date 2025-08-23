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
    
    # Contribution to log partition function (not free energy)
    # According to polymer theory: log(Z) = log(Z_BP) + Σ φ(W) Z_W
    contribution = phi_W * Z_W 
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
η = 0.4
nsamples = 1
w = 8

cluster_data = load_latest_cluster_file(N,w)
clusters_by_site = cluster_data.clusters_by_site
all_loops = cluster_data.all_loops

tensors, peps = peps_controllable(N, T; η=η, ti=ti, orthog=orthog)
exact_FE_density = (log(real(contract_peps_no_phys(peps; cutoff=1E-12, maxdim=2^N)))) / (N*T)

adj_mat, edges, links = BP.get_adj_mat(tensors)
messages = BP.get_messages(tensors, edges, links; random_part=0.0)
messages = BP.message_passing(tensors, messages, edges, adj_mat; α = annealing, max_iters=maxiter)
println("self consistency chec...", BP.check_self_consistency(tensors, messages, adj_mat))
Z_l = BP.get_fixed_point_list(tensors, messages, adj_mat)
bp_FE_density = sum(log.(real.(Z_l)))  / (N*T)
T_normalized = BP.normalize_tensors(tensors, Z_l)


global clustercorrx = 0

# Option 2: Deduplicate clusters globally using Set to automatically remove duplicates
println("Collecting all unique clusters from all sites...")
all_unique_clusters = Set{Cluster}()
global total_clusters_before_dedup = 0

# Debug: Let's also collect clusters manually to analyze them
clusters_by_signature = Dict{Tuple, Vector{Tuple{Int, Cluster}}}()

for site in 1:(N*T)
    clusters = clusters_by_site[site]
    global total_clusters_before_dedup += length(clusters)
    for cluster in clusters
        push!(all_unique_clusters, cluster)  # Set automatically deduplicates
        
        # Debug: Create a signature for manual comparison
        signature = (cluster.weight, cluster.total_loops, sort(cluster.loop_ids), sort(collect(cluster.multiplicities)))
        if !haskey(clusters_by_signature, signature)
            clusters_by_signature[signature] = []
        end
        push!(clusters_by_signature[signature], (site, cluster))
    end
end

println("Before deduplication: $total_clusters_before_dedup total cluster instances")
println("After Set deduplication: $(length(all_unique_clusters)) unique clusters")
println("Manual signature analysis: $(length(clusters_by_signature)) unique signatures")

# Let's see if there are actually any duplicates based on our manual signatures
global manual_duplicates = 0
for (sig, cluster_list) in clusters_by_signature
    if length(cluster_list) > 1
        global manual_duplicates += length(cluster_list) - 1
    end
end
println("Manual analysis found $manual_duplicates duplicates that should be removed")

# Compute contributions for each unique cluster exactly once
println("Computing contributions...")
if manual_duplicates > 0
    println("Using manual deduplication (Set didn't work properly)")
    # Use manual deduplication - take one representative from each signature group
    for (sig, cluster_list) in clusters_by_signature
        cluster = cluster_list[1][2]  # Take the first cluster from this signature group
        contribution = cluster_contr(T_normalized, messages, edges, links, adj_mat, cluster, all_loops)
        # Only add contribution if computation was successful (non-zero)
        if !isnan(contribution) && isfinite(contribution)
            global clustercorrx += contribution
        end
    end
else
    println("Using Set deduplication (seems to be working)")
    for cluster in all_unique_clusters
        contribution = cluster_contr(T_normalized, messages, edges, links, adj_mat, cluster, all_loops)
        # Only add contribution if computation was successful (non-zero)
        if !isnan(contribution) && isfinite(contribution)
            global clustercorrx += contribution
        end
    end
end

# The correction is now the sum over unique clusters (no multiplication needed)
cluster_FE_density_correction = clustercorrx / (N*T) 

println("BP error: ", abs(bp_FE_density - exact_FE_density))
println("BP + cluster correction error: ", abs(bp_FE_density + cluster_FE_density_correction - exact_FE_density))