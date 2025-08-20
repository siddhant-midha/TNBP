#!/usr/bin/env julia

"""
test_cluster_expansion_debug_new.jl

Debug the cluster expansion calculation to understand why cluster corrections 
are oscillating with large values. Focus on individual contributions.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("../functions/Ising2D.jl")
include("../functions/BP.jl")
using Serialization
using ITensors, Graphs

function debug_cluster_expansion_single_beta()
    """
    Debug cluster expansion for a single β value to understand the large corrections.
    """
    
    println("🔍 Debugging Cluster Expansion Calculation")
    println("="^60)
    
    # Parameters
    L = 11
    β = 0.5  # Choose a moderate β where we see large corrections
    h = 0.0
    
    println("Parameters: L=$L, β=$β, h=$h")
    println()
    
    # Step 1: Create Ising tensor network
    println("Step 1: Creating Ising tensor network...")
    T = Ising2D.get_ising_tn(L, β; h=h)
    N = L^2
    println("✅ Created $N tensors for $(L)×$(L) lattice")
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    println("✅ Found $(length(edges)) edges")
    
    # Step 2: Compute BP fixed point
    println("\nStep 2: Computing BP fixed point...")
    messages = BP.get_messages(T, edges, links)
    messages = BP.message_passing(T, messages, edges, adj_mat; max_iters=1000)
    println("✅ BP converged")
    
    # Get BP partition function (before normalization)
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Z_bp_full)
    println("✅ BP log partition function: $log_Z_bp_full")
    
    # Step 3: Normalize tensors
    println("\nStep 3: Normalizing tensors...")
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)
    
    # Step 4: Load saved clusters
    cluster_file = "../saved_clusters/periodic_clusters_L11_w10_2025-08-14T20-25-54-769.jld2"
    println("\nStep 4: Loading saved clusters...")
    
    loaded_data = open(cluster_file, "r") do io
        deserialize(io)
    end
    
    data = loaded_data["data"]
    summary = loaded_data["summary"]
    
    # Get clusters from saved data
    if hasfield(typeof(data), :clusters_by_site)
        all_clusters = data.clusters_by_site[1]  # Site 1 clusters only
    else
        all_clusters = data.clusters
    end
    
    all_loops = data.all_loops
    
    println("✅ Loaded $(length(all_clusters)) clusters from site 1")
    
    # Step 5: Debug weight-4 cluster contributions
    println("\nStep 5: Debugging weight-4 cluster contributions...")
    weight4_clusters = [c for c in all_clusters if c.weight == 4]
    
    println("Found $(length(weight4_clusters)) weight-4 clusters on site 1")
    
    total_correction = 0.0 + 0im
    
    for (i, cluster) in enumerate(weight4_clusters)
        println("\n--- Cluster $i ---")
        println("Total loops: $(cluster.total_loops)")
        println("Loop IDs: $(cluster.loop_ids)")
        println("Weight: $(cluster.weight)")
        
        # Compute Ursell function φ(W)
        phi_W = 1.0  # For single loop clusters, φ(W) = 1
        println("Ursell function φ(W): $phi_W")
        
        # Compute cluster correction Z_W = Z_l (since single loop with multiplicity 1)
        Z_W = 1.0 + 0im
        
        for loop_id in cluster.loop_ids
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            
            println("Loop $loop_id:")
            println("  Multiplicity: $multiplicity")
            println("  Loop weight: $(loop.weight)")
            println("  Loop vertices: $(loop.vertices)")
            println("  Loop edges: $(loop.edges)")
            
            # Convert loop edges format for BP.jl
            loop_edges_bp = Tuple{Int,Int}[]
            for edge in loop.edges
                v1, v2 = edge
                push!(loop_edges_bp, (min(v1, v2), max(v1, v2)))
            end
            println("  BP edges: $loop_edges_bp")
            
            # Compute loop contribution Z_l
            try
                Z_l_tensor = BP.loop_contribution(loop_edges_bp, messages, T_normalized, edges, links, adj_mat)
                println("  Z_l tensor: $Z_l_tensor")
                println("  Z_l tensor indices: $(inds(Z_l_tensor))")
                
                if isempty(inds(Z_l_tensor))
                    Z_l_scalar = scalar(Z_l_tensor)
                    println("  Z_l scalar: $Z_l_scalar")
                    println("  |Z_l|: $(abs(Z_l_scalar))")
                    
                    Z_W *= Z_l_scalar^multiplicity
                    println("  Z_W after this loop: $Z_W")
                else
                    println("  ⚠️  Z_l tensor is not scalar! Has indices: $(inds(Z_l_tensor))")
                end
                
            catch e
                println("  ❌ Error computing loop contribution: $e")
                continue
            end
        end
        
        # Add weighted contribution: φ(W) * Z_W
        contribution = phi_W * Z_W
        total_correction += contribution
        
        println("Cluster contribution: φ(W) * Z_W = $phi_W * $Z_W = $contribution")
        println("Running total correction: $total_correction")
        
        # Check if contribution is reasonable
        if abs(contribution) > 1.0
            println("⚠️  Large contribution detected!")
        end
        if abs(imag(contribution)) > 1e-10
            println("⚠️  Non-negligible imaginary part!")
        end
    end
    
    println("\n" * "="^60)
    println("Final Results:")
    println("Total correction from $(length(weight4_clusters)) weight-4 clusters: $total_correction")
    println("Real part: $(real(total_correction))")
    println("Imaginary part: $(imag(total_correction))")
    
    # Compare with what we expect
    println("\nAnalysis:")
    if abs(total_correction) > 10.0
        println("❌ Correction is suspiciously large")
    elseif abs(total_correction) > 1.0
        println("⚠️  Correction is moderately large")
    else
        println("✅ Correction magnitude seems reasonable")
    end
    
    # Also compute what the corrected free energy would be
    corrected_log_Z = log_Z_bp_full + total_correction
    bp_free_energy = -real(log_Z_bp_full) / (2 * N)
    corrected_free_energy = -real(log_Z_bp_full) / (2 * N) + real(total_correction)
    
    println("\nFree Energy Analysis:")
    println("BP free energy density: $bp_free_energy")
    println("Corrected free energy density: $corrected_free_energy")
    println("Correction to free energy density: $(real(total_correction))")
    
    # Compare with exact
    exact_free_energy = Ising2D.free_energy(β)
    println("Exact free energy (Onsager): $exact_free_energy")
    println("BP error: $(abs(bp_free_energy - exact_free_energy))")
    println("Corrected error: $(abs(corrected_free_energy - exact_free_energy))")
    
    return total_correction
end

# Run the debug
if abspath(PROGRAM_FILE) == @__FILE__
    debug_cluster_expansion_single_beta()
end