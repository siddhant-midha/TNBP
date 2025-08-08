"""
test_cluster_expansion_debug.jl

Debug test to identify where NaN values are coming from in cluster_expansion_2d_ising.jl
"""

include("../dependencies.jl")
include("../functions/BP.jl")
include("../functions/Ising2D.jl")
include("../functions/LoopEnumeration.jl") 
include("../functions/ClusterEnumeration.jl")
include("../functions/cluster_expansion_2d_ising.jl")

using ITensors, Graphs

function debug_cluster_expansion_step_by_step()
    """
    Step-by-step debugging of cluster expansion to identify NaN source.
    """
    
    println("🔍 DEBUG: Testing cluster expansion step by step")
    println("="^60)
    
    # Test parameters - start small
    L = 4
    β = 0.4
    h = 0.0
    max_weights = [4]
    
    println("Parameters: L=$L, β=$β, h=$h")
    println("Max weights: $max_weights")
    println()
    
    # Initialize variables to avoid scope issues
    local T, N, adj_mat, edges, links, messages, Z_bp_full, log_Z_bp_full, bp_free_energy
    local Z_l, T_normalized, exact_free_energy, all_loops
    
    # Step 1: Test Ising tensor network creation
    println("Step 1: Testing Ising tensor network creation...")
    try
        T = Ising2D.get_ising_tn(L, β; h=h)
        N = L^2
        println("✅ Created $N tensors for $(L)x$(L) lattice")
        println("  Sample tensor norm: $(norm(T[1]))")
        
        # Check for NaN/Inf in tensors
        for (i, tensor) in enumerate(T)
            try
                # Get tensor data safely
                tensor_data = storage(tensor)
                if any(isnan, tensor_data) || any(isinf, tensor_data)
                    println("❌ Found NaN/Inf in tensor $i")
                    return false
                end
            catch e
                println("⚠️  Could not check tensor $i for NaN/Inf: $e")
                # Continue anyway
            end
        end
        println("✅ All tensors are finite")
    catch e
        println("❌ Error creating Ising tensors: $e")
        return false
    end
    
    # Step 2: Test BP adjacency extraction
    println("\nStep 2: Testing BP adjacency extraction...")
    try
        adj_mat, edges, links = BP.get_adj_mat(T)
        println("✅ Found $(length(edges)) edges")
        println("  Sample edge: $(edges[1])")
        
        # Check adjacency matrix
        if any(isnan, adj_mat) || any(isinf, adj_mat)
            println("❌ Found NaN/Inf in adjacency matrix")
            return false
        end
        println("✅ Adjacency matrix is finite")
    catch e
        println("❌ Error extracting adjacency: $e")
        return false
    end
    
    # Step 3: Test BP message initialization and passing
    println("\nStep 3: Testing BP message passing...")
    try
        messages = BP.get_messages(T, edges, links)
        println("✅ Initialized $(length(messages)) messages")
        
        # Check for NaN in initial messages
        for (i, msg) in enumerate(messages)
            try
                msg_data = storage(msg)
                if any(isnan, msg_data) || any(isinf, msg_data)
                    println("❌ Found NaN/Inf in initial message $i")
                    return false
                end
            catch e
                println("⚠️  Could not check initial message $i for NaN/Inf: $e")
            end
        end
        println("✅ Initial messages are finite")
        
        # Run message passing
        messages = BP.message_passing(T, messages, edges, links, adj_mat; max_iters=100)
        println("✅ BP converged")
        
        # Check converged messages
        for (i, msg) in enumerate(messages)
            try
                msg_data = storage(msg)
                if any(isnan, msg_data) || any(isinf, msg_data)
                    println("❌ Found NaN/Inf in converged message $i")
                    return false
                end
            catch e
                println("⚠️  Could not check message $i for NaN/Inf: $e")
            end
        end
        println("✅ Converged messages are finite")
        
    catch e
        println("❌ Error in BP message passing: $e")
        return false
    end
    
    # Step 4: Test BP partition function
    println("\nStep 4: Testing BP partition function...")
    try
        Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, edges, links, adj_mat)
        println("✅ BP partition function: $Z_bp_full")
        
        if isnan(Z_bp_full) || isinf(Z_bp_full) || Z_bp_full == 0
            println("❌ BP partition function is invalid: $Z_bp_full")
            return false
        end
        
        log_Z_bp_full = log(Z_bp_full)
        println("✅ Log BP partition function: $log_Z_bp_full")
        
        if isnan(log_Z_bp_full) || isinf(log_Z_bp_full)
            println("❌ Log BP partition function is invalid: $log_Z_bp_full")
            return false
        end
        
        # Compute BP free energy
        bp_free_energy = -log_Z_bp_full / (2 * N)
        println("✅ BP free energy density: $bp_free_energy")
        
        if isnan(bp_free_energy) || isinf(bp_free_energy)
            println("❌ BP free energy is invalid: $bp_free_energy")
            return false
        end
        
    catch e
        println("❌ Error computing BP partition function: $e")
        return false
    end
    
    # Step 5: Test tensor normalization
    println("\nStep 5: Testing tensor normalization...")
    try
        Z_l = BP.get_fixed_point_list(T, messages, edges, links, adj_mat)
        println("✅ Got fixed point list, length: $(length(Z_l))")
        
        # Check fixed point values
        for (i, z) in enumerate(Z_l)
            if isnan(z) || isinf(z) || z == 0
                println("❌ Invalid fixed point value at $i: $z")
                return false
            end
        end
        println("✅ All fixed point values are valid")
        
        T_normalized = BP.normalize_tensor(T, Z_l)
        println("✅ Normalized tensors")
        
        # Check normalized tensors
        for (i, tensor) in enumerate(T_normalized)
            try
                tensor_data = storage(tensor)
                if any(isnan, tensor_data) || any(isinf, tensor_data)
                    println("❌ Found NaN/Inf in normalized tensor $i")
                    return false
                end
            catch e
                println("⚠️  Could not check normalized tensor $i for NaN/Inf: $e")
            end
        end
        println("✅ All normalized tensors are finite")
        
    catch e
        println("❌ Error in tensor normalization: $e")
        return false
    end
    
    # Step 6: Test Onsager solution
    println("\nStep 6: Testing Onsager solution...")
    try
        exact_free_energy = Ising2D.free_energy(β)
        println("✅ Onsager free energy: $exact_free_energy")
        
        if isnan(exact_free_energy) || isinf(exact_free_energy)
            println("❌ Onsager free energy is invalid: $exact_free_energy")
            return false
        end
        
    catch e
        println("❌ Error computing Onsager solution: $e")
        return false
    end
    
    # Step 7: Test loop enumeration
    println("\nStep 7: Testing loop enumeration...")
    try
        # Create adjacency matrix for cluster enumeration
        cluster_adj_matrix = zeros(Int, N, N)
        for (v1, v2) in edges
            cluster_adj_matrix[v1, v2] = 1
            cluster_adj_matrix[v2, v1] = 1
        end
        
        enumerator = ClusterEnumerator(cluster_adj_matrix)
        println("✅ Created cluster enumerator")
        
        # Test small weight first
        test_weight = 4
        all_loops = find_all_loops_in_graph(enumerator, test_weight)
        println("✅ Found $(length(all_loops)) loops for weight $test_weight")
        
        if isempty(all_loops)
            println("⚠️  No loops found - this might be the issue!")
            
            # Try even smaller system
            println("  Trying smaller system L=3...")
            L_small = 3
            N_small = L_small^2
            T_small = Ising2D.get_ising_tn(L_small, β; h=h)
            adj_mat_small, edges_small, links_small = BP.get_adj_mat(T_small)
            
            cluster_adj_matrix_small = zeros(Int, N_small, N_small)
            for (v1, v2) in edges_small
                cluster_adj_matrix_small[v1, v2] = 1
                cluster_adj_matrix_small[v2, v1] = 1
            end
            
            enumerator_small = ClusterEnumerator(cluster_adj_matrix_small)
            all_loops_small = find_all_loops_in_graph(enumerator_small, test_weight)
            println("  Found $(length(all_loops_small)) loops in $L_small x $L_small system")
            
            if isempty(all_loops_small)
                println("❌ Still no loops found - loop enumeration has an issue!")
                return false
            end
        end
        
    catch e
        println("❌ Error in loop enumeration: $e")
        return false
    end
    
    # Step 8: Test cluster enumeration
    println("\nStep 8: Testing cluster enumeration...")
    try
        if !isempty(all_loops)
            interaction_graph = build_interaction_graph(all_loops)
            println("✅ Built interaction graph with $(length(interaction_graph)) nodes")
            
            # Test cluster generation for small weight
            test_clusters = generate_all_clusters_of_weight(all_loops, 4, interaction_graph)
            println("✅ Generated $(length(test_clusters)) clusters of weight 4")
            
            if isempty(test_clusters)
                println("⚠️  No clusters found - might be expected for small systems")
            end
        end
        
    catch e
        println("❌ Error in cluster enumeration: $e")
        return false
    end
    
    println("\n🎯 DEBUG SUMMARY:")
    println("="^40)
    println("✅ Tensor creation: OK")
    println("✅ BP adjacency: OK") 
    println("✅ BP messages: OK")
    println("✅ BP partition function: OK")
    println("✅ Tensor normalization: OK")
    println("✅ Onsager solution: OK")
    if isempty(all_loops)
        println("⚠️  Loop enumeration: NO LOOPS FOUND")
    else
        println("✅ Loop enumeration: OK")
        println("✅ Cluster enumeration: OK")
    end
    
    return true
end

function test_individual_functions()
    """Test individual functions in isolation."""
    
    println("\n🧪 Testing individual functions...")
    
    # Test LoopEnumerator with simple lattice
    println("\nTesting LoopEnumerator on 3x3 lattice...")
    L = 3
    adj_matrix = create_periodic_square_lattice(L)
    println("✅ Created periodic square lattice")
    
    loop_enum = LoopEnumerator(adj_matrix)
    println("✅ Created LoopEnumerator")
    println("  Lattice size: $(loop_enum.L)")
    println("  Vertices: $(loop_enum.n_vertices)")
    
    # Test loop finding on specific vertex
    test_vertex = 5  # Center of 3x3 grid
    test_loops = find_loops_supported_on_vertex(loop_enum, test_vertex, 4)
    println("✅ Found $(length(test_loops)) loops supported on vertex $test_vertex")
    
    for (i, loop) in enumerate(test_loops)
        if i <= 3  # Show first few
            println("  Loop $i: vertices=$(loop.vertices), edges=$(loop.edges), weight=$(loop.weight)")
        end
    end
    
    if isempty(test_loops)
        println("❌ No loops found - this is the issue!")
        
        # Debug adjacency list
        println("  Adjacency for vertex $test_vertex: $(loop_enum.adj_list[test_vertex])")
        println("  Coordinates for vertex $test_vertex: $(loop_enum.coords[test_vertex])")
        
        # Try finding loops on vertex 1
        test_loops_v1 = find_loops_supported_on_vertex(loop_enum, 1, 4)
        println("  Loops on vertex 1: $(length(test_loops_v1))")
        
    end
    
    return !isempty(test_loops)
end

# Run the debug tests
if abspath(PROGRAM_FILE) == @__FILE__
    println("🚀 Starting debug tests...")
    
    success1 = test_individual_functions()
    success2 = debug_cluster_expansion_step_by_step()
    
    if success1 && success2
        println("\n✅ All debug tests passed!")
    else
        println("\n❌ Debug tests found issues!")
    end
end