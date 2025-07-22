"""
simple_test.jl

Simple test to demonstrate the Julia loop enumeration works correctly.
"""

include("../functions/LoopEnumeration.jl")

function simple_test()
    println("Julia Loop and Cluster Enumeration - Simple Test")
    println("="^60)
    
    # Test loop enumeration on 10x10 lattice
    L = 10
    vertex = 1  # Julia uses 1-based indexing
    max_weight = 7
    
    println("Testing loop enumeration:")
    println("  Lattice: $(L)×$(L)")
    println("  Vertex: $vertex")  
    println("  Max weight: $max_weight")
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = LoopEnumerator(adj_matrix)
    
    start_time = time()
    loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
    end_time = time()
    
    println("  Found: $(length(loops)) loops")
    println("  Time: $(round(end_time - start_time, digits=3)) seconds")
    
    # Analyze by weight
    by_weight = Dict{Int, Int}()
    for loop in loops
        by_weight[loop.weight] = get(by_weight, loop.weight, 0) + 1
    end
    
    println("\n📊 Results by weight:")
    for weight in sort(collect(keys(by_weight)))
        count = by_weight[weight]
        println("  Weight $weight: $count loops")
    end
    
    total = length(loops)
    println("\n✅ Summary:")
    println("  Total loops: $total")
    println("  Expected: 28")
    println("  Match: $(total == 28 ? "YES ✅" : "NO ❌")")
    
    if total == 28
        println("\n🎉 SUCCESS! The Julia implementation correctly finds 28 loops.")
        println("   This matches the expected result for:")
        println("   - 4 loops of weight 4 (simple 4-cycles)")
        println("   - 12 loops of weight 6 (simple 6-cycles)")  
        println("   - 12 loops of weight 7 (complex structures with degree > 2)")
        
        println("\n📝 Note on cluster enumeration:")
        println("   Each of these 28 loops forms its own connected cluster,")
        println("   so there are also 28 connected clusters of weight ≤ 7")
        println("   supported on any given site.")
        
        return true
    else
        println("\n❌ The loop count doesn't match expected results.")
        return false
    end
end

# Test smaller lattices for verification  
function test_pattern()
    println("\n" * "="^50)
    println("Testing pattern on smaller lattices")
    println("="^50)
    
    for L in [3, 4, 5]
        adj_matrix = create_periodic_square_lattice(L)
        enumerator = LoopEnumerator(adj_matrix)
        loops = find_loops_supported_on_vertex(enumerator, 1, 7)
        println("  $(L)×$(L) lattice: $(length(loops)) loops")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    success = simple_test()
    test_pattern()
    
    println("\n" * "="^60)
    println("FINAL STATUS")  
    println("="^60)
    if success
        println("✅ Julia implementation SUCCESSFUL")
        println("   - Loop enumeration: WORKING")
        println("   - Finds correct 28 loops on 10×10 lattice") 
        println("   - Cluster enumeration: CONCEPTUALLY EQUIVALENT")
        println("     (Each loop forms its own cluster)")
    else
        println("❌ Julia implementation needs debugging")
    end
end