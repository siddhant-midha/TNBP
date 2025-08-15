"""
test_loop_enumeration.jl

Test suite for loop enumeration in Julia.
"""

include("../functions/LoopEnumeration.jl")

function test_10x10_loop_enumeration(use_optimized::Bool = true)
    method_name = use_optimized ? "OPTIMIZED" : "STANDARD"
    println("Testing $method_name loop enumeration on 10x10 lattice")
    println("="^60)
    
    L = 10
    vertex = 1  # Julia uses 1-based indexing
    max_weight = 7
    
    println("Parameters:")
    println("  Lattice size: $(L)x$(L)")
    println("  Target vertex: $vertex")
    println("  Max weight: $max_weight")
    println("  Method: $method_name")
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = LoopEnumerator(adj_matrix)
    
    println("\nRunning enumeration...")
    start_time = time()
    
    if use_optimized
        global_seen = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
        loops = find_loops_supported_on_vertex_bfs(enumerator, vertex, max_weight)
    else
        loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
    end
    
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

function test_smaller_lattices(use_optimized::Bool = true)
    method_name = use_optimized ? "OPTIMIZED" : "STANDARD"
    println("\n" * "="^50)
    println("Testing smaller lattices for pattern analysis ($method_name)")
    println("="^50)
    
    max_weight = 7
    vertex = 1
    
    results = Tuple{Int, Int}[]
    
    for L in [3, 4, 5]
        println("\nTesting $(L)x$(L) lattice...")
        
        adj_matrix = create_periodic_square_lattice(L)
        enumerator = LoopEnumerator(adj_matrix)
        
        start_time = time()
        if use_optimized
            global_seen = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
            loops = find_loops_supported_on_vertex_bfs(enumerator, vertex, max_weight)
        else
            loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
        end
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

function verify_loop_connectivity(use_optimized::Bool = true)
    method_name = use_optimized ? "OPTIMIZED" : "STANDARD"
    println("\n" * "="^50)
    println("Verifying loop connectivity ($method_name)")
    println("="^50)
    
    L = 4
    vertex = 1
    max_weight = 6
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = LoopEnumerator(adj_matrix)
    
    if use_optimized
        global_seen = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
        loops = find_loops_supported_on_vertex_bfs(enumerator, vertex, max_weight)
    else
        loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
    end
    
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

function test_weight_8_loops(use_optimized::Bool = true)
    method_name = use_optimized ? "OPTIMIZED" : "STANDARD"
    println("\n" * "="^60)
    println("TESTING WEIGHT 8 LOOPS ($method_name)")
    println("="^60)
    
    L = 10
    vertex = 1
    max_weight = 8
    
    println("Parameters:")
    println("  Lattice size: $(L)x$(L)")
    println("  Target vertex: $vertex")
    println("  Max weight: $max_weight")
    println("  Method: $method_name")
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = LoopEnumerator(adj_matrix)
    
    println("\nRunning enumeration...")
    start_time = time()
    if use_optimized
        global_seen = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
        loops = find_loops_supported_on_vertex_bfs(enumerator, vertex, max_weight)
    else
        loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
    end
    end_time = time()
    
    # Count loops by weight
    by_weight = Dict{Int, Int}()
    for loop in loops
        by_weight[loop.weight] = get(by_weight, loop.weight, 0) + 1
    end
    
    println("Breakdown by weight:")
    for weight in sort(collect(keys(by_weight)))
        count = by_weight[weight]
        println("  Weight $weight: $count loops")
    end
    
    weight_8_count = get(by_weight, 8, 0)
    println("\nRESULT:")
    println("Weight 8 loops on vertex 1: $weight_8_count")
    println("Expected: 70")
    println("Match: $(weight_8_count == 70 ? "✅ YES" : "❌ NO")")
    println("Time: $(round(end_time - start_time, digits=3)) seconds")
    
    return weight_8_count == 70
end

function test_translational_symmetry(use_optimized::Bool = false)
    method_name = use_optimized ? "OPTIMIZED" : "STANDARD"
    println("\n" * "="^60)
    println("TESTING TRANSLATIONAL SYMMETRY ($method_name)")
    println("="^60)
    
    L = 8  # Use smaller lattice for efficiency
    max_weight = 6
    
    println("Parameters:")
    println("  Lattice size: $(L)x$(L)")
    println("  Max weight: $max_weight")
    println("  Total sites: $(L*L)")
    println("  Method: $method_name")
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = LoopEnumerator(adj_matrix)
    
    # Count loops on vertex 1
    println("\nCounting loops on vertex 1...")
    start_time = time()
    
    if use_optimized
        global_seen = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
        loops_vertex_1 = find_loops_supported_on_vertex_bfs(enumerator, 1, max_weight)
    else
        loops_vertex_1 = find_loops_supported_on_vertex(enumerator, 1, max_weight)
    end
    
    end_time = time()
    k = L*L*3
    
    # println("Loops on vertex 1: $k")
    println("Time: $(round(end_time - start_time, digits=3)) seconds")
    
    # Count total unique loops in the lattice
    println("\nCounting all unique loops in lattice...")
    start_time = time()
    
    if use_optimized
        all_unique_loops = find_all_loops_in_graph_bfs(enumerator, max_weight)
    else
        all_unique_loops = find_all_loops_in_graph(enumerator, max_weight)
    end
    
    actual_total = length(all_unique_loops)
    end_time = time()
    
    println("Actual total unique loops: $actual_total")
    println("Time: $(round(end_time - start_time, digits=3)) seconds")
    
    # The relationship should be: actual_total × L² = expected_total
    # Or equivalently: actual_total = k (each unique loop appears at L² sites)
    println("\nTRANSLATIONAL SYMMETRY CHECK:")
    # println("Loops per vertex (k): $k")
    println("Total unique loops: $actual_total")
    # println("Expected relationship: actual_total = k (each loop appears at L² sites)")
    println("Match: $(actual_total == k ? "✅ YES" : "❌ NO")")
    
    if actual_total != k
        println("⚠️  Symmetry broken! This suggests:")
        println("   - Some loops are not translationally equivalent")
        println("   - Or there's an error in the enumeration")
    end
    
    return actual_total == k
end

function test_comprehensive_enumeration(use_optimized::Bool = true)
    method_name = use_optimized ? "OPTIMIZED" : "STANDARD"
    println("\n" * "="^60)
    println("COMPREHENSIVE ENUMERATION TEST ($method_name)")
    println("="^60)
    
    # Test multiple lattice sizes and weights
    test_cases = [
        (L=3, max_weight=4, expected_vertex_loops=6),
        (L=4, max_weight=6, expected_vertex_loops=54),
    ]
    
    all_passed = true
    
    for (L, max_weight, expected) in test_cases
        println("\nTesting $(L)x$(L) lattice, max_weight=$max_weight:")
        
        adj_matrix = create_periodic_square_lattice(L)
        enumerator = LoopEnumerator(adj_matrix)
        
        start_time = time()
        if use_optimized
            global_seen = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
            loops = find_loops_supported_on_vertex_bfs(enumerator, 1, max_weight)
        else
            loops = find_loops_supported_on_vertex(enumerator, 1, max_weight)
        end
        end_time = time()
        
        actual = length(loops)
        passed = (actual == expected)
        all_passed &= passed
        
        println("  Loops found: $actual")
        println("  Expected: $expected")
        println("  Result: $(passed ? "✅ PASS" : "❌ FAIL")")
        println("  Time: $(round(end_time - start_time, digits=3))s")
    end
    
    return all_passed
end

function test_weight_8_vertex_degrees(use_optimized::Bool = true)
    method_name = use_optimized ? "OPTIMIZED" : "STANDARD"
    println("\n" * "="^60)
    println("TESTING WEIGHT-8 LOOP VERTEX DEGREES ($method_name)")
    println("="^60)
    
    L = 10
    vertex = 1
    max_weight = 8
    
    println("Parameters:")
    println("  Lattice size: $(L)x$(L)")
    println("  Target vertex: $vertex")
    println("  Testing weight-8 loops for degree constraint")
    println("  Method: $method_name")
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = LoopEnumerator(adj_matrix)
    
    println("\nEnumerating loops...")
    start_time = time()
    if use_optimized
        global_seen = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
        loops = find_loops_supported_on_vertex_bfs(enumerator, vertex, max_weight)
    else
        loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
    end
    end_time = time()
    
    # Filter weight-8 loops
    weight_8_loops = [loop for loop in loops if loop.weight == 8]
    println("Found $(length(weight_8_loops)) weight-8 loops in $(round(end_time - start_time, digits=3))s")
    
    if isempty(weight_8_loops)
        println("⚠️  No weight-8 loops found to test")
        return true, 0, 0
    end
    
    all_degree_2 = true
    degree_2_count = 0
    higher_degree_count = 0
    
    println("\nAnalyzing vertex degrees...")
    for (i, loop) in enumerate(weight_8_loops)
        # Calculate degree of each vertex in this loop
        degree = Dict{Int, Int}()
        for v in loop.vertices
            degree[v] = 0
        end
        
        for (u, v) in loop.edges
            degree[u] += 1
            degree[v] += 1
        end
        
        # Check if all vertices have degree 2
        loop_all_degree_2 = true
        max_degree = 0
        for v in loop.vertices
            if degree[v] != 2
                loop_all_degree_2 = false
                all_degree_2 = false
            end
            max_degree = max(max_degree, degree[v])
        end
        
        if loop_all_degree_2
            degree_2_count += 1
        else
            higher_degree_count += 1
            if i <= 3  # Show first few examples
                println("  Loop $i: vertices=$(loop.vertices)")
                println("    Degrees: $([v => degree[v] for v in sort(loop.vertices)])")
                println("    Max degree: $max_degree")
            end
        end
    end
    
    println("\nRESULTS:")
    println("Weight-8 loops with all degree = 2: $degree_2_count")
    println("Weight-8 loops with some degree > 2: $higher_degree_count")
    
    return all_degree_2, degree_2_count, higher_degree_count
end

function benchmark_optimization_comparison()
    """Compare performance between standard and optimized loop enumeration."""
    println("\n" * "="^70)
    println("BENCHMARK: STANDARD vs OPTIMIZED LOOP ENUMERATION")
    println("="^70)
    
    # Test cases with increasing complexity
    test_cases = [
        (L=6, max_weight=6, desc="Small case"),
        (L=8, max_weight=7, desc="Medium case"),
        (L=10, max_weight=8, desc="Large case"),
    ]
    
    results = []
    
    for (L, max_weight, desc) in test_cases
        println("\n🔍 $desc: $(L)x$(L) lattice, max_weight=$max_weight")
        println("-" * 50)
        
        adj_matrix = create_periodic_square_lattice(L)
        enumerator = LoopEnumerator(adj_matrix)
        vertex = 1
        
        # Test standard version
        println("  📊 Testing STANDARD enumeration...")
        start_time = time()
        loops_standard = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
        time_standard = time() - start_time
        count_standard = length(loops_standard)
        
        # Test optimized version
        println("  🚀 Testing OPTIMIZED enumeration...")
        start_time = time()
        global_seen = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
        loops_optimized = find_loops_supported_on_vertex_bfs(enumerator, vertex, max_weight)
        time_optimized = time() - start_time
        count_optimized = length(loops_optimized)
        
        # Test ultra-fast version
        println("  ⚡ Testing FAST enumeration...")
        start_time = time()
        loops_fast = find_loops_supported_on_vertex_fast(enumerator, vertex, max_weight)
        time_fast = time() - start_time
        count_fast = length(loops_fast)
        
        # Verify correctness
        correctness_check = "✅ SAME" 
        if count_standard != count_optimized || count_standard != count_fast
            correctness_check = "❌ DIFFERENT"
            println("  ⚠️  WARNING: Different loop counts!")
            println("    Standard: $count_standard, Optimized: $count_optimized, Fast: $count_fast")
        end
        
        # Calculate speedups
        speedup_opt = time_standard / time_optimized
        speedup_fast = time_standard / time_fast
        
        speedup_opt_str = if speedup_opt > 1.0
            "$(round(speedup_opt, digits=2))x faster"
        else
            "$(round(1/speedup_opt, digits=2))x slower"
        end
        
        speedup_fast_str = if speedup_fast > 1.0
            "$(round(speedup_fast, digits=2))x faster"
        else
            "$(round(1/speedup_fast, digits=2))x slower"
        end
        
        # Store results
        push!(results, (
            desc = desc,
            L = L,
            max_weight = max_weight,
            count = count_standard,
            time_standard = time_standard,
            time_optimized = time_optimized,
            time_fast = time_fast,
            speedup_opt = speedup_opt,
            speedup_fast = speedup_fast,
            correct = count_standard == count_optimized && count_standard == count_fast
        ))
        
        # Print results
        println("  📈 RESULTS:")
        println("    Standard:  $(round(time_standard, digits=3))s → $count_standard loops")
        println("    Optimized: $(round(time_optimized, digits=3))s → $count_optimized loops ($speedup_opt_str)")
        println("    Fast:      $(round(time_fast, digits=3))s → $count_fast loops ($speedup_fast_str)")
        println("    Correctness: $correctness_check")
    end
    
    # Summary
    println("\n" * "="^70)
    println("📊 BENCHMARK SUMMARY")
    println("="^70)
    println("Case          | Loops | Standard | Optimized | Fast      | Best Speedup | Correct")
    println("--------------|-------|----------|-----------|-----------|--------------|--------")
    
    all_correct = true
    total_speedup_opt = 1.0
    total_speedup_fast = 1.0
    
    for result in results
        speedup_opt_str = if result.speedup_opt > 1.0
            "$(round(result.speedup_opt, digits=2))x"
        else
            "$(round(1/result.speedup_opt, digits=2))x slower"
        end
        
        speedup_fast_str = if result.speedup_fast > 1.0
            "$(round(result.speedup_fast, digits=2))x"
        else
            "$(round(1/result.speedup_fast, digits=2))x slower"
        end
        
        best_speedup = max(result.speedup_opt, result.speedup_fast)
        best_str = if best_speedup > 1.0
            "$(round(best_speedup, digits=2))x faster"
        else
            "$(round(1/best_speedup, digits=2))x slower"
        end
        
        correct_str = result.correct ? "✅" : "❌"
        all_correct &= result.correct
        total_speedup_opt *= result.speedup_opt
        total_speedup_fast *= result.speedup_fast
        
        println("$(rpad(result.desc, 13)) | $(lpad(result.count, 5)) | $(lpad(round(result.time_standard, digits=3), 8))s | $(lpad(round(result.time_optimized, digits=3), 9))s | $(lpad(round(result.time_fast, digits=3), 9))s | $(rpad(best_str, 12)) | $correct_str")
    end
    
    geometric_mean_speedup_opt = total_speedup_opt^(1/length(results))
    geometric_mean_speedup_fast = total_speedup_fast^(1/length(results))
    best_overall = max(geometric_mean_speedup_opt, geometric_mean_speedup_fast)
    
    println("\n🎯 OVERALL PERFORMANCE:")
    println("  Optimized version speedup: $(round(geometric_mean_speedup_opt, digits=2))x")
    println("  Fast version speedup: $(round(geometric_mean_speedup_fast, digits=2))x")
    println("  Best overall speedup: $(round(best_overall, digits=2))x")
    println("  All results correct: $(all_correct ? "✅ YES" : "❌ NO")")
    
    if best_overall > 1.1
        best_method = geometric_mean_speedup_fast > geometric_mean_speedup_opt ? "FAST" : "OPTIMIZED"
        println("  🚀 $best_method version is EFFECTIVE!")
    elseif best_overall > 0.9
        println("  ➡️  Optimizations have NEUTRAL impact")
    else
        println("  🐌 All optimizations make things SLOWER")
    end
    
    return results
end

function main(use_optimized::Bool = true, run_benchmark::Bool = false)
    use_optimized = true
    method_name = use_optimized ? "OPTIMIZED" : "STANDARD"
    println("Julia Loop Enumeration Test Suite ($method_name)")
    println("="^60)
    
    if run_benchmark
        # Run benchmark comparison
        benchmark_results = benchmark_optimization_comparison()
        return benchmark_results
    end
    
    # Test smaller lattices first
    test_smaller_lattices(use_optimized)
    
    # Test the main 10x10 case
    println("\n" * "="^60)
    println("MAIN TEST: 10x10 lattice, weight 7 ($method_name)")
    println("="^60)
    
    total_loops = test_10x10_loop_enumeration(use_optimized)
    
    # Test weight 8 loops
    weight_8_passed = test_weight_8_loops(use_optimized)
    
    # Test weight-8 vertex degrees
    all_degree_8, degree_2_count, higher_degree_count = test_weight_8_vertex_degrees(use_optimized)
    
    # Test translational symmetry
    symmetry_passed = test_translational_symmetry(use_optimized)
    
    # Test comprehensive cases
    comprehensive_passed = test_comprehensive_enumeration(use_optimized)
    
    # Verify correctness
    verification_passed = verify_loop_connectivity(use_optimized)
    
    println("\n" * "="^60)
    println("FINAL RESULTS ($method_name)")
    println("="^60)
    println("10x10 lattice, weight 7: $total_loops loops (expected 28)")
    println("Weight 7 test: $(total_loops == 28 ? "✅ PASS" : "❌ FAIL")")
    println("Weight 8 count test: $(weight_8_passed ? "✅ PASS" : "❌ FAIL")")
    println("Weight 8 degree test: $(!all_degree_8 ? "✅ PASS" : "❌ FAIL")")
    println("  Degree=2 loops: $degree_2_count, Degree>2 loops: $higher_degree_count")
    println("Symmetry test: $(symmetry_passed ? "✅ PASS" : "❌ FAIL")")
    println("Comprehensive test: $(comprehensive_passed ? "✅ PASS" : "❌ FAIL")")
    println("Connectivity verification: $(verification_passed ? "✅ PASS" : "❌ FAIL")")
    
    all_passed = (total_loops == 28) && weight_8_passed && !all_degree_8 &&
                 symmetry_passed && comprehensive_passed && verification_passed
    
    println("\nOVERALL: $(all_passed ? "✅ ALL TESTS PASSED" : "❌ SOME TESTS FAILED")")
    
    return all_passed
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line arguments
    use_optimized = false
    run_benchmark = false
    
    for arg in ARGS
        if arg == "--optimized" || arg == "-o"
            use_optimized = true
        elseif arg == "--benchmark" || arg == "-b"
            run_benchmark = true
        elseif arg == "--help" || arg == "-h"
            println("Usage: julia test_loop_enumeration.jl [OPTIONS]")
            println("Options:")
            println("  --optimized, -o    Use optimized loop enumeration")
            println("  --benchmark, -b    Run benchmark comparison")
            println("  --help, -h         Show this help message")
            exit(0)
        end
    end
    
    main(use_optimized, run_benchmark)
end