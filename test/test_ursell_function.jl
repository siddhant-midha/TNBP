"""
test_ursell_function.jl

Test the Ursell function implementation for both connected and disconnected clusters.
"""

include("../functions/ClusterEnumeration.jl")

function main()
    println("="^60)
    println("Testing Ursell Function Implementation")
    println("="^60)
    
    # Run the basic tests from ClusterEnumeration.jl
    test_ursell_function()
    
    println("\n" * "="^60)
    println("Additional Ursell Function Tests")
    println("="^60)
    
    # Test with a smaller lattice for more predictable results
    L = 3
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = ClusterEnumerator(adj_matrix)
    
    # Find all loops up to weight 6
    all_loops = find_all_loops_in_graph(enumerator, 6)
    println("Found $(length(all_loops)) loops on $(L)×$(L) lattice")
    
    # Test various cluster configurations
    if length(all_loops) >= 1
        println("\n📋 Testing different cluster types:")
        
        # Test 1: Single loop cluster (should give φ = 1)
        loop1 = all_loops[1]
        single_cluster = Cluster([1], Dict(1 => 1), loop1.weight, 1)
        ursell_val = ursell_function(single_cluster, [loop1])
        println("  Single loop cluster: φ = $ursell_val (expected: 1.0)")
        
        # Test 2: Double loop cluster with same loop (multiplicity 2)
        if loop1.weight * 2 <= 10  # Keep total weight reasonable
            double_cluster = Cluster([1], Dict(1 => 2), loop1.weight * 2, 2)
            ursell_val = ursell_function(double_cluster, [loop1])
            println("  Same loop twice (η=2): φ = $ursell_val (expected: 1/2! = 0.5)")
        end
        
        # Test 3: Two different loops (if available)
        if length(all_loops) >= 2
            loop2 = all_loops[2]
            
            # Check if loops are compatible (don't share vertices)
            are_compatible = !loops_are_incompatible(loop1, loop2)
            
            two_loop_cluster = Cluster([1, 2], Dict(1 => 1, 2 => 1), 
                                     loop1.weight + loop2.weight, 2)
            ursell_val = ursell_function(two_loop_cluster, [loop1, loop2])
            
            if are_compatible
                println("  Two compatible loops: φ = $ursell_val (expected: 0.0 - disconnected)")
            else
                println("  Two incompatible loops: φ = $ursell_val (expected: non-zero - connected)")
            end
        end
    end
    
    println("\n" * "="^60)
    println("Mathematical Properties Test")
    println("="^60)
    
    # Verify key mathematical properties
    test_mathematical_properties(all_loops)
    
    println("\n✅ Ursell function testing completed!")
end

function test_mathematical_properties(all_loops::Vector{Loop})
    """Test key mathematical properties of the Ursell function."""
    
    if length(all_loops) < 2
        println("⚠️  Need at least 2 loops for mathematical property tests")
        return
    end
    
    println("Testing mathematical properties:")
    
    # Property 1: Single loop always gives φ = 1
    loop1 = all_loops[1]
    single_cluster = Cluster([1], Dict(1 => 1), loop1.weight, 1)
    phi_single = ursell_function(single_cluster, [loop1])
    
    if abs(phi_single - 1.0) < 1e-10
        println("  ✅ Property 1: Single loop φ = 1")
    else
        println("  ❌ Property 1: Single loop φ = $phi_single ≠ 1")
    end
    
    # Property 2: Disconnected clusters give φ = 0
    found_disconnected = false
    for i in 1:min(3, length(all_loops))
        for j in (i+1):min(4, length(all_loops))
            loop_i = all_loops[i]
            loop_j = all_loops[j]
            
            if !loops_are_incompatible(loop_i, loop_j)
                # Found disconnected pair
                disconnected_cluster = Cluster([i, j], Dict(i => 1, j => 1), 
                                             loop_i.weight + loop_j.weight, 2)
                phi_disc = ursell_function(disconnected_cluster, all_loops)
                
                if abs(phi_disc) < 1e-10
                    println("  ✅ Property 2: Disconnected cluster φ = 0")
                else
                    println("  ❌ Property 2: Disconnected cluster φ = $phi_disc ≠ 0")
                end
                found_disconnected = true
                break
            end
        end
        if found_disconnected
            break
        end
    end
    
    if !found_disconnected
        println("  ⚠️  Property 2: No disconnected loops found for testing")
    end
    
    # Property 3: Multiplicity affects the result (φ scales with 1/η!)
    if loop1.weight * 2 <= 8  # Reasonable total weight
        eta1_cluster = Cluster([1], Dict(1 => 1), loop1.weight, 1)
        eta2_cluster = Cluster([1], Dict(1 => 2), loop1.weight * 2, 2)
        
        phi_eta1 = ursell_function(eta1_cluster, [loop1])
        phi_eta2 = ursell_function(eta2_cluster, [loop1])
        
        expected_ratio = 1.0 / factorial(2)  # 1/2!
        actual_ratio = phi_eta2 / phi_eta1
        
        if abs(actual_ratio - expected_ratio) < 1e-10
            println("  ✅ Property 3: Multiplicity scaling φ(η=2)/φ(η=1) = 1/2!")
        else
            println("  ❌ Property 3: Expected ratio 1/2! = $expected_ratio, got $actual_ratio")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end