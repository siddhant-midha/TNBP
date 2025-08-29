include("../../functions/ClusterEnumeration.jl")
include("../../functions/Ising2D.jl")
include("../../functions/BP.jl")
include("../../test/test_utils.jl")
include("../../functions/brute_force_new.jl")

using Serialization
using Plots
using ITensors, Graphs

# Data structure for single-site enumeration results (copied from generate_ising_clusters_one_site.jl)
struct SingleSiteClusterData
    """Data structure to save single-site cluster enumeration results."""
    clusters::Vector{Cluster}
    all_loops::Vector{Loop}
    adj_matrix::Matrix{Int}
    max_weight::Int
    lattice_size::Int
    site::Int
    enumeration_time::Float64
    translation_removal_time::Float64
    timestamp::String
    boundary_condition::String
    clusters_before_translation_removal::Int
    clusters_after_translation_removal::Int
    canonical_forms_count::Int
end


function cluster_expansion_2d_ising_with_single_site_clusters(L, max_weight, T_normalized, messages, edges, links, adj_mat, filename)
    println("📖 Loading cluster data from: $(filename)")
    loaded_data = open(filename, "r") do io
        deserialize(io)
    end
    cluster_data = loaded_data["data"]

    cluster_contributions = compute_all_single_site_cluster_contributions(
        T_normalized, messages, edges, links, adj_mat, cluster_data, L, max_weight)
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

function compute_all_single_site_cluster_contributions(T_normalized, messages, edges, links, adj_mat, 
                                                     cluster_data, L::Int, max_weight::Int)
    """
    Compute ALL cluster contributions once for single-site cluster data.
    Returns a vector of dictionaries with contribution details for each cluster.
    """
    
    # Get clusters from single-site data
    all_clusters = cluster_data.clusters
    all_loops = cluster_data.all_loops
    
    println("📊 Computing contributions for all single-site clusters")
    println("   Available clusters: $(length(all_clusters))")
    println("   Available loops: $(length(all_loops))")
    
    # Filter clusters by weight
    relevant_clusters = [c for c in all_clusters if c.weight <= max_weight]
    println("   Using $(length(relevant_clusters)) clusters (weight ≤ $max_weight)")
    
    if isempty(relevant_clusters)
        println("⚠️  No relevant clusters found")
        return []
    end
    
    # Compute contributions for all relevant clusters
    cluster_contributions = []
    
    println("💫 Computing individual cluster contributions...")
    
    for (i, cluster) in enumerate(relevant_clusters)
        if i % 100 == 0 || i <= 10
            println("  Processing cluster $i/$(length(relevant_clusters)) (weight=$(cluster.weight))...")
        end
        
        # Compute Ursell function φ(W)
        phi_W = ursell_function(cluster, all_loops)
        
        if abs(phi_W) < 1e-15
            continue  # Skip negligible contributions
        end
        
        # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
        Z_W = 1.0 + 0im  # Complex number for cluster contribution
        computation_successful = true
        
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
                println("⚠️  Error computing loop contribution for loop $loop_id: $e")
                computation_successful = false
                break
            end
        end
        
        if !computation_successful
            continue
        end
        
        # The contribution to log(Z) is φ(W) * Z_W, so contribution to f is -φ(W) * Z_W / 2.0
        contribution = -phi_W * Z_W #/ 2.0  # Factor of 2 comes from Ising model convention -> Not sure about this
        
        # Store contribution details
        contrib_info = Dict(
            "cluster_id" => i,
            "weight" => cluster.weight,
            "phi_W" => phi_W,
            "Z_W" => Z_W,
            "contribution" => real(contribution)
        )
        push!(cluster_contributions, contrib_info)
        
        if abs(contribution) > 1e-10 && i <= 20
            println("    Cluster $i (weight=$(cluster.weight)): φ=$phi_W, Z_W=$Z_W, f contribution=$contribution")
        end
    end
    
    println("✅ Computed $(length(cluster_contributions)) valid cluster contributions")
    
    # Sort by weight for easier analysis
    sort!(cluster_contributions, by=x->x["weight"])
    
    return cluster_contributions
end
