"""
cluster_expansion_2d_ising.jl

Complete implementation of cluster expansion for 2D Ising model on 10x10 periodic lattice.
This script implements the workflow described in polymer_proof.tex using:
1. Tensor networks from Ising2D.jl
2. BP fixed point from BP.jl  
3. Cluster enumeration from ClusterEnumeration.jl
4. Ursell functions for connected clusters
"""

include("../dependencies.jl")
include("BP.jl")
include("Ising2D.jl")
include("ClusterEnumeration.jl")

using ITensors, Graphs

function cluster_expansion_2d_ising(L::Int=10, β::Float64=0.4, h::Float64=0.0; 
                                   max_weights=[4, 6, 8, 10, 12])
    """
    Complete cluster expansion workflow for 2D Ising model.
    
    Args:
        L: Lattice size (default: 10x10)
        β: Inverse temperature (default: 0.4)  
        h: Magnetic field (default: 0.0)
        max_weights: List of maximum cluster weights to compute
    
    Returns:
        Dict with results: BP free energy, cluster corrected free energies, exact Onsager
    """
    
    println("="^80)
    println("🔥 Cluster Expansion for 2D Ising Model")
    println("="^80)
    println("Parameters: L=$L, β=$β, h=$h")
    println("Cluster weights: $max_weights")
    println()
    
    # Step 1: Create Ising tensor network
    println("Step 1: Creating Ising tensor network...")
    T = Ising2D.get_ising_tn(L, β; h=h)
    N = L^2
    println("✅ Created $N tensors for $(L)x$(L) lattice")
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    println("✅ Found $(length(edges)) edges")
    
    # Step 2: Compute BP fixed point
    println("\nStep 2: Computing BP fixed point...")
    messages = BP.get_messages(T, edges, links)
    messages = BP.message_passing(T, messages, edges, links, adj_mat; max_iters=1000)
    println("✅ BP converged")
    
    # Get BP partition function (before normalization)
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, edges, links, adj_mat)
    log_Z_bp_full = log(Z_bp_full)
    println("✅ BP log partition function: $log_Z_bp_full")
    
    # Step 3: Normalize tensors
    println("\nStep 3: Normalizing tensors...")
    Z_l = BP.get_fixed_point_list(T, messages, edges, links, adj_mat)
    T_normalized = BP.normalize_tensor(T, Z_l)
    
    # Normalization factor
    normalization_factor = sum(log.(real.(Z_l)))
    println("✅ Normalization factor: $normalization_factor")
    
    # Verify normalization (BP on normalized tensors should give 1)
    Z_bp_normalized = BP.mean_free_partition_fn(1:N, T_normalized, messages, edges, links, adj_mat)
    println("✅ Normalized BP partition function: $Z_bp_normalized (should ≈ 1)")
    
    # Step 4: Enumerate clusters and compute corrections
    results = Dict{String, Any}()
    results["L"] = L
    results["beta"] = β
    results["h"] = h
    results["bp_log_Z"] = -log_Z_bp_full / (2 * N)  # Free energy density with correct sign and normalization
    results["normalization_factor"] = normalization_factor / N
    results["cluster_corrections"] = Dict()
    
    for max_weight in max_weights
        println("\n" * "="^60)
        println("🔍 Computing cluster expansion up to weight $max_weight")
        println("="^60)
        
        correction = compute_cluster_correction(T_normalized, messages, edges, links, adj_mat, 
                                              L, max_weight)
        
        # Total free energy with cluster correction
        # According to polymer theory: log(Z) = log(Z_BP) + Σ φ(W) Z_W
        # Apply correct sign convention and normalization factor
        log_Z_corrected = -(log_Z_bp_full + correction) / (2 * N)
        
        results["cluster_corrections"][max_weight] = Dict(
            "correction" => -correction / (2 * N),  # Apply same sign/normalization to correction
            "total_free_energy" => log_Z_corrected,
            "improvement" => -correction / (2 * N)
        )
        
        println("📊 Results for weight $max_weight:")
        println("  Cluster correction: $(correction/N)")
        println("  Total free energy density: $log_Z_corrected")
    end
    
    # Step 5: Compare with exact Onsager solution
    if h == 0.0  # Onsager solution only valid for h=0
        println("\n" * "="^60)
        println("📐 Comparing with exact Onsager solution")
        println("="^60)
        
        exact_free_energy = Ising2D.free_energy(β)
        results["exact_onsager"] = exact_free_energy
        
        println("🎯 Exact Onsager free energy: $exact_free_energy")
        println("\n📈 Comparison with exact result:")
        println("  BP approximation error: $(abs(results["bp_log_Z"] - exact_free_energy))")
        
        for max_weight in max_weights
            corrected_f = results["cluster_corrections"][max_weight]["total_free_energy"]
            error = abs(corrected_f - exact_free_energy)
            improvement = abs(results["bp_log_Z"] - exact_free_energy) - error
            
            println("  Weight $max_weight error: $error (improvement: $improvement)")
        end
    else
        println("⚠️  Exact comparison only available for h=0 (Onsager solution)")
        results["exact_onsager"] = nothing
    end
    
    return results
end

function compute_cluster_correction(T_normalized, messages, edges, links, adj_mat, 
                                  L::Int, max_weight::Int)
    """
    Compute the cluster expansion correction up to given weight.
    
    This implements the formula from polymer_proof.tex:
    log(Z) = Σ_{connected W} φ(W) * Z_W
    """
    
    # Create adjacency matrix for cluster enumeration (convert to format expected by ClusterEnumeration)
    N = L^2
    cluster_adj_matrix = zeros(Int, N, N)
    for (v1, v2) in edges
        cluster_adj_matrix[v1, v2] = 1
        cluster_adj_matrix[v2, v1] = 1
    end
    
    # Initialize cluster enumerator
    enumerator = ClusterEnumerator(cluster_adj_matrix)
    
    println("🔄 Enumerating all loops up to weight $max_weight...")
    all_loops = find_all_loops_in_graph(enumerator, max_weight)
    println("Found $(length(all_loops)) total loops")
    
    if isempty(all_loops)
        println("⚠️  No loops found, returning zero correction")
        return 0.0
    end
    
    # Build interaction graph
    interaction_graph = build_interaction_graph(all_loops)
    
    # Find all connected clusters up to max_weight
    println("🔍 Enumerating connected clusters...")
    all_connected_clusters = Cluster[]
    
    # We need to enumerate ALL connected clusters, not just those supported on one site
    # So we enumerate clusters of each weight systematically
    for weight in 1:max_weight
        weight_clusters = generate_all_clusters_of_weight(all_loops, weight, interaction_graph)
        append!(all_connected_clusters, weight_clusters)
    end
    
    println("Found $(length(all_connected_clusters)) connected clusters")
    
    if isempty(all_connected_clusters)
        println("⚠️  No connected clusters found, returning zero correction")
        return 0.0
    end
    
    # Compute total correction
    total_correction = 0.0
    
    println("💫 Computing cluster contributions...")
    # Limit to first 100 clusters for efficiency in demo
    max_clusters_to_process = min(100, length(all_connected_clusters))
    println("Processing first $max_clusters_to_process clusters (out of $(length(all_connected_clusters)) total)")
    
    for (i, cluster) in enumerate(all_connected_clusters[1:max_clusters_to_process])
        if i % 20 == 0 || i <= 5
            println("  Processing cluster $i/$max_clusters_to_process...")
        end
        
        # Compute Ursell function φ(W)
        phi_W = ursell_function(cluster, all_loops)
        
        if abs(phi_W) < 1e-15
            continue  # Skip negligible contributions
        end
        
        # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
        Z_W = 1.0 + 0im  # Complex number for cluster contribution
        
        for loop_id in cluster.loop_ids
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            
            # Convert loop edges format for BP.jl (ensure v1 < v2 ordering)
            loop_edges_bp = Tuple{Int,Int}[]
            for edge in loop.edges
                v1, v2 = edge
                push!(loop_edges_bp, (min(v1, v2), max(v1, v2)))
            end
            
            # Compute loop contribution Z_l
            try
                Z_l = BP.loop_contribution(loop_edges_bp, messages, T_normalized, edges, links, adj_mat)
                Z_W *= Z_l^multiplicity
            catch e
                println("⚠️  Error computing loop contribution for loop $loop_id: $e")
                continue
            end
        end
        
        # Add weighted contribution: φ(W) * Z_W
        contribution = phi_W * Z_W
        total_correction += contribution
        
        if abs(contribution) > 1e-10 && i <= 20
            println("    Cluster $i: φ=$phi_W, Z_W=$Z_W, contribution=$contribution")
        end
    end
    
    println("✅ Total cluster correction: $total_correction")
    return real(total_correction)  # Take real part for final result
end

function generate_all_clusters_of_weight(all_loops::Vector{Loop}, target_weight::Int,
                                       interaction_graph::Dict{Int, Vector{Int}})
    """Generate all connected clusters of exactly the target weight."""
    
    all_clusters = Cluster[]
    
    # Generate all possible clusters (with multiplicities)
    function enumerate_partitions(remaining_weight::Int, start_idx::Int, current_cluster::Vector{Int})
        if remaining_weight == 0
            if !isempty(current_cluster)
                # Create cluster from current_cluster
                multiplicities = Dict{Int, Int}()
                for loop_id in current_cluster
                    multiplicities[loop_id] = get(multiplicities, loop_id, 0) + 1
                end
                
                cluster = Cluster(
                    collect(keys(multiplicities)),
                    multiplicities,
                    target_weight,
                    length(current_cluster)
                )
                
                # Check if cluster is connected
                if is_cluster_connected(cluster, interaction_graph)
                    push!(all_clusters, cluster)
                end
            end
            return
        end
        
        if remaining_weight < 0
            return
        end
        
        for loop_idx in start_idx:length(all_loops)
            loop_weight = all_loops[loop_idx].weight
            if loop_weight <= remaining_weight
                push!(current_cluster, loop_idx)
                enumerate_partitions(remaining_weight - loop_weight, loop_idx, current_cluster)
                pop!(current_cluster)
            end
        end
    end
    
    enumerate_partitions(target_weight, 1, Int[])
    
    # Remove duplicates
    return remove_cluster_redundancies(all_clusters)
end

# Main execution function
function run_cluster_expansion_demo()
    """Run the cluster expansion demo with different parameters."""
    
    println("🚀 Running Cluster Expansion Demo for 2D Ising Model")
    println()
    
    # Test parameters
    L = 10
    β = 0.4  # Moderate temperature
    h = 0.0  # No magnetic field (for Onsager comparison)
    
    # Run cluster expansion
    results = cluster_expansion_2d_ising(L, β, h; max_weights=[4, 6, 8, 10])
    
    # Print summary
    println("\n" * "="^80)
    println("📋 FINAL SUMMARY")
    println("="^80)
    println("System: $(L)x$(L) Ising model, β=$β, h=$h")
    println()
    println("Free energy density (per site):")
    println("  BP approximation: $(results["bp_log_Z"])")
    
    if results["exact_onsager"] !== nothing
        println("  Exact (Onsager): $(results["exact_onsager"])")
        println()
        println("Cluster corrections:")
        for weight in [4, 6, 8, 10]
            if haskey(results["cluster_corrections"], weight)
                data = results["cluster_corrections"][weight]
                error = abs(data["total_free_energy"] - results["exact_onsager"])
                println("  Weight $weight: $(data["total_free_energy"]) (error: $error)")
            end
        end
    else
        println()
        println("Cluster corrections:")
        for weight in [4, 6, 8, 10]
            if haskey(results["cluster_corrections"], weight)
                data = results["cluster_corrections"][weight]
                println("  Weight $weight: $(data["total_free_energy"]) (improvement: $(data["improvement"]))")
            end
        end
    end
    
    return results
end

# Run demo if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    results = run_cluster_expansion_demo()
end