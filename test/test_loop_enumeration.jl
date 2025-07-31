"""
test_loop_enumeration.jl

Test suite for loop enumeration in Julia.
"""

include("../functions/LoopEnumeration.jl")

function test_10x10_loop_enumeration()
    println("Testing loop enumeration on 10x10 lattice")
    println("="^60)
    
    L = 10
    vertex = 1  # Julia uses 1-based indexing
    max_weight = 7
    
    println("Parameters:")
    println("  Lattice size: $(L)x$(L)")
    println("  Target vertex: $vertex")
    println("  Max weight: $max_weight")
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = LoopEnumerator(adj_matrix)
    
    println("\nRunning enumeration...")
    start_time = time()
    loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
    end_time = time()
    
    total_loops = length(loops)
    println("Found $total_loops loops in $(round(end_time - start_time, digits=3)) seconds")
    
    # Analyze by weight
    by_weight = Dict{Int, Vector{Loop}}()
    for loop in loops
        if !haskey(by_weight, loop.weight)
            by_weight[loop.weight] = Loop[]
        end
        push!(by_weight[loop.weight], loop)
    end
    
    println("\nBreakdown by weight:")
    println("Weight | Count | Examples")
    println("-------|-------|----------")
    
    for weight in sort(collect(keys(by_weight)))
        loops_of_weight = by_weight[weight]
        count = length(loops_of_weight)
        
        # Show examples
        examples = String[]
        for i in 1:min(2, count)
            loop = loops_of_weight[i]
            vertices = loop.vertices
            edges = loop.edges
            
            # Count degree of each vertex in this loop
            degree = Dict{Int, Int}()
            for v in vertices
                degree[v] = 0
            end
            for (u, v) in edges
                degree[u] += 1
                degree[v] += 1
            end
            
            # Find vertices with degree > 2
            high_degree = [v for v in vertices if degree[v] > 2]
            
            if !isempty(high_degree)
                push!(examples, "V:$(length(vertices)), deg>2:$high_degree")
            else
                push!(examples, "V:$(length(vertices)), all deg=2")
            end
        end
        
        examples_str = join(examples, "; ")
        if count > 2
            examples_str *= "; +$(count-2) more"
        end
        
        println("  $(lpad(weight, 2))   | $(lpad(count, 4))  | $examples_str")
    end
    
    println("\nSUMMARY:")
    println("Total loops with weight <= $max_weight: $total_loops")
    println("Expected: 28")
    println("Match: $(total_loops == 28 ? "✅ YES" : "❌ NO")")
    
    return total_loops
end

function test_smaller_lattices()
    println("\n" * "="^50)
    println("Testing smaller lattices for pattern analysis")
    println("="^50)
    
    max_weight = 7
    vertex = 1
    
    results = Tuple{Int, Int}[]
    
    for L in [3, 4, 5]
        println("\nTesting $(L)x$(L) lattice...")
        
        adj_matrix = create_periodic_square_lattice(L)
        enumerator = LoopEnumerator(adj_matrix)
        
        start_time = time()
        loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
        end_time = time()
        
        total = length(loops)
        push!(results, (L, total))
        
        # Quick breakdown
        by_weight = Dict{Int, Int}()
        for loop in loops
            by_weight[loop.weight] = get(by_weight, loop.weight, 0) + 1
        end
        
        weight_str = join(["w$w:$c" for (w, c) in sort(collect(by_weight))], ", ")
        println("  $(L)x$(L): $total total ($weight_str) in $(round(end_time - start_time, digits=3))s")
    end
    
    println("\nPattern analysis:")
    for (L, total) in results
        println("  $(L)x$(L) lattice: $total loops")
    end
end

function verify_loop_connectivity()
    println("\n" * "="^50)
    println("Verifying loop connectivity")
    println("="^50)
    
    L = 4
    vertex = 1
    max_weight = 6
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = LoopEnumerator(adj_matrix)
    
    loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
    
    println("Testing connectivity of found loops...")
    all_connected = true
    
    for (i, loop) in enumerate(loops)
        vertices_set = Set(loop.vertices)
        if !is_connected(vertices_set, loop.edges)
            println("❌ Loop $i is not connected!")
            all_connected = false
        end
    end
    
    if all_connected
        println("✅ All $(length(loops)) loops are properly connected")
    end
    
    # Test degree constraint
    println("\nTesting degree >= 2 constraint...")
    all_valid_degree = true
    
    for (i, loop) in enumerate(loops)
        vertices_set = Set(loop.vertices)
        if !is_valid_loop(vertices_set, loop.edges)
            println("❌ Loop $i violates degree >= 2 constraint!")
            all_valid_degree = false
        end
    end
    
    if all_valid_degree
        println("✅ All $(length(loops)) loops satisfy degree >= 2 constraint")
    end
    
    return all_connected && all_valid_degree
end

function main()
    println("Julia Loop Enumeration Test Suite")
    println("="^60)
    
    # Test smaller lattices first
    test_smaller_lattices()
    
    # Test the main 10x10 case
    println("\n" * "="^60)
    println("MAIN TEST: 10x10 lattice")
    println("="^60)
    
    total_loops = test_10x10_loop_enumeration()
    
    # Verify correctness
    verification_passed = verify_loop_connectivity()
    
    println("\n" * "="^60)
    println("FINAL RESULT")
    println("="^60)
    println("10x10 lattice, max weight 7: $total_loops loops")
    println("Expected: 28 loops")
    println("Match: $(total_loops == 28 ? "✅ YES" : "❌ NO")")
    println("Verification: $(verification_passed ? "✅ PASSED" : "❌ FAILED")")
    
    return total_loops == 28 && verification_passed
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end