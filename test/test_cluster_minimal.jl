"""
test_cluster_minimal.jl

Test cluster expansion directly to find the NaN source
"""

include("../dependencies.jl")
include("../functions/BP.jl")
include("../functions/Ising2D.jl")
include("../functions/LoopEnumeration.jl") 
include("../functions/ClusterEnumeration.jl")
include("../functions/cluster_expansion_2d_ising.jl")

using ITensors, Graphs

function test_cluster_expansion_minimal()
    """Test cluster expansion on minimal system."""
    
    println("🧪 Testing cluster expansion...")
    
    L = 4
    β = 0.4
    h = 0.0
    max_weights = [4]
    
    try
        println("Calling cluster_expansion_2d_ising...")
        results = cluster_expansion_2d_ising(L, β, h; max_weights=max_weights)
        
        println("Results keys: $(keys(results))")
        println("BP free energy: $(results["bp_log_Z"])")
        println("Exact free energy: $(results["exact_onsager"])")
        
        if haskey(results, "cluster_corrections")
            println("Cluster corrections: $(keys(results["cluster_corrections"]))")
            for weight in max_weights
                if haskey(results["cluster_corrections"], weight)
                    corr_data = results["cluster_corrections"][weight]
                    println("  Weight $weight: $(corr_data)")
                else
                    println("  Weight $weight: MISSING")
                end
            end
        else
            println("No cluster_corrections in results!")
        end
        
        return true
        
    catch e
        println("❌ Error in cluster expansion: $e")
        println("Error type: $(typeof(e))")
        if isa(e, BoundsError)
            println("This is a bounds error - likely indexing issue")
        elseif isa(e, UndefRefError)
            println("This is an undefined reference error")
        elseif isa(e, MethodError)
            println("This is a method error - likely function signature issue")
        end
        
        # Print stack trace for more info
        println("Stack trace:")
        for (i, frame) in enumerate(stacktrace(catch_backtrace()))
            if i <= 10  # Show first 10 frames
                println("  $i: $frame")
            end
        end
        
        return false
    end
end

function test_cluster_functions_individually()
    """Test cluster functions step by step."""
    
    println("\n🧪 Testing cluster functions individually...")
    
    L = 3  # Even smaller
    β = 0.4
    h = 0.0
    
    try
        # Step 1: Create system
        T = Ising2D.get_ising_tn(L, β; h=h)
        N = L^2
        adj_mat, edges, links = BP.get_adj_mat(T)
        messages = BP.get_messages(T, edges, links)
        messages = BP.message_passing(T, messages, edges, links, adj_mat; max_iters=100)
        println("✅ BP setup complete")
        
        # Step 2: Test cluster enumeration setup
        cluster_adj_matrix = zeros(Int, N, N)
        for (v1, v2) in edges
            cluster_adj_matrix[v1, v2] = 1
            cluster_adj_matrix[v2, v1] = 1
        end
        
        enumerator = ClusterEnumerator(cluster_adj_matrix)
        println("✅ Created cluster enumerator")
        
        # Step 3: Test loop finding
        all_loops = find_all_loops_in_graph(enumerator, 4)
        println("✅ Found $(length(all_loops)) loops")
        
        if isempty(all_loops)
            println("⚠️  No loops found - this could be the issue!")
            return false
        end
        
        # Step 4: Test interaction graph
        interaction_graph = build_interaction_graph(all_loops)
        println("✅ Built interaction graph")
        
        # Step 5: Test cluster generation
        test_clusters = generate_all_clusters_of_weight(all_loops, 4, interaction_graph)
        println("✅ Generated $(length(test_clusters)) clusters")
        
        if isempty(test_clusters)
            println("⚠️  No clusters found")
            return false
        end
        
        # Step 6: Test normalization
        Z_l = BP.get_fixed_point_list(T, messages, edges, links, adj_mat)
        T_normalized = BP.normalize_tensor(T, Z_l)
        println("✅ Normalized tensors")
        
        # Step 7: Test loop contribution calculation
        if !isempty(test_clusters)
            cluster = test_clusters[1]
            println("Testing cluster: $(cluster)")
            
            for loop_id in cluster.loop_ids
                loop = all_loops[loop_id]
                loop_edges_bp = [(min(e[1], e[2]), max(e[1], e[2])) for e in loop.edges]
                
                try
                    Z_l = BP.loop_contribution(loop_edges_bp, messages, T_normalized, edges, links, adj_mat)
                    println("✅ Loop $loop_id contribution: $Z_l")
                    
                    if isnan(Z_l) || isinf(Z_l)
                        println("❌ Loop contribution is invalid: $Z_l")
                        return false
                    end
                    
                catch e
                    println("❌ Error computing loop contribution: $e")
                    return false
                end
            end
        end
        
        return true
        
    catch e
        println("❌ Error in individual testing: $e")
        return false
    end
end

# Run tests
if abspath(PROGRAM_FILE) == @__FILE__
    println("🚀 Testing cluster expansion in detail...")
    
    success1 = test_cluster_functions_individually()
    success2 = test_cluster_expansion_minimal()
    
    if success1 && success2
        println("\n✅ All cluster expansion tests passed")
    else
        println("\n❌ Found issues in cluster expansion")
    end
end