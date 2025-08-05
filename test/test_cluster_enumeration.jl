"""
test_cluster_enumeration.jl

Test suite for cluster enumeration in Julia.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

function test_10x10_weight_7_cluster_enumeration()
    println("Testing cluster enumeration on 10x10 lattice")
    println("(Allowing vertices with degree > 2)")
    println("="^60)
    
    L = 10
    site = 1  # Julia uses 1-based indexing
    max_weight = 7
    
    println("Parameters:")
    println("  Lattice size: $(L)x$(L)")
    println("  Target site: $site")
    println("  Max cluster weight: $max_weight")
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = ClusterEnumerator(adj_matrix)
    
    println("\nStarting enumeration...")
    start_time = time()
    
    try
        clusters = enumerate_connected_clusters(enumerator, site, max_weight)
        end_time = time()
        
        println("\n✅ SUCCESS!")
        println("Found $(length(clusters)) connected clusters in $(round(end_time - start_time, digits=3)) seconds")
        
        analyze_cluster_results(clusters, max_weight)
        return length(clusters)
        
    catch e
        end_time = time()
        println("\n❌ ERROR after $(round(end_time - start_time, digits=3)) seconds:")
        println("   $e")
        return nothing
    end
end


function test_10x10_weight_8_cluster_enumeration()
    println("Testing cluster enumeration on 10x10 lattice")
    println("(Allowing vertices with degree > 2)")
    println("="^60)
    
    L = 10
    site = 1  # Julia uses 1-based indexing
    max_weight = 8
    
    println("Parameters:")
    println("  Lattice size: $(L)x$(L)")
    println("  Target site: $site")
    println("  Max cluster weight: $max_weight")
    
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = ClusterEnumerator(adj_matrix)
    
    println("\nStarting enumeration...")
    start_time = time()
    
    try
        clusters = enumerate_connected_clusters(enumerator, site, max_weight)
        end_time = time()
        
        println("\n✅ SUCCESS!")
        println("Found $(length(clusters)) connected clusters in $(round(end_time - start_time, digits=3)) seconds")
        
        analyze_cluster_results(clusters, max_weight)
        return length(clusters)
        
    catch e
        end_time = time()
        println("\n❌ ERROR after $(round(end_time - start_time, digits=3)) seconds:")
        println("   $e")
        return nothing
    end
end


function analyze_cluster_results(clusters::Vector{Cluster}, max_weight::Int)
    """Analyze and display cluster enumeration results."""
    
    # Group by total weight
    by_weight = Dict{Int, Vector{Cluster}}()
    for cluster in clusters
        if !haskey(by_weight, cluster.weight)
            by_weight[cluster.weight] = Cluster[]
        end
        push!(by_weight[cluster.weight], cluster)
    end
    
    println("\n📊 Analysis by cluster weight:")
    println("Weight | Count | Examples")
    println("-------|-------|----------")
    
    for weight in sort(collect(keys(by_weight)))
        clusters_of_weight = by_weight[weight]
        count = length(clusters_of_weight)
        
        # Create examples
        examples = String[]
        for i in 1:min(2, count)
            cluster = clusters_of_weight[i]
            
            # Show multiplicities
            mult_items = String[]
            total_loop_instances = 0
            for loop_id in cluster.loop_ids
                mult = cluster.multiplicities[loop_id]
                push!(mult_items, "L$(loop_id)×$(mult)")
                total_loop_instances += mult
            end
            
            mult_str = join(mult_items, ",")
            push!(examples, "($mult_str)")
        end
        
        examples_str = join(examples, "; ")
        if count > 2
            examples_str *= "; +$(count-2) more"
        end
        
        println("  $(lpad(weight, 2))   | $(lpad(count, 4))  | $examples_str")
    end
    
    # Analyze cluster composition
    println("\n🔍 Cluster composition analysis:")
    
    # Count single-loop vs multi-loop clusters
    single_loop_clusters = sum(1 for c in clusters if length(c.loop_ids) == 1)
    multi_loop_clusters = length(clusters) - single_loop_clusters
    
    println("  Single-loop clusters: $single_loop_clusters")
    println("  Multi-loop clusters: $multi_loop_clusters")
    
    # Analyze multiplicities
    max_multiplicity = 0
    multiplicity_counts = Dict{Int, Int}()
    
    for cluster in clusters
        for (loop_id, mult) in cluster.multiplicities
            max_multiplicity = max(max_multiplicity, mult)
            multiplicity_counts[mult] = get(multiplicity_counts, mult, 0) + 1
        end
    end
    
    println("  Maximum loop multiplicity: $max_multiplicity")
    print("  Multiplicity distribution: ")
    mult_dist = join(["$(mult)×$(count)" for (mult, count) in sort(collect(multiplicity_counts))], ", ")
    println(mult_dist)
    
    # Show some detailed examples
    println("\n📝 Detailed examples:")
    
    example_count = 0
    for weight in sort(collect(keys(by_weight)))
        if example_count >= 5  # Show at most 5 examples
            break
        end
        
        clusters_of_weight = by_weight[weight]
        
        for cluster in clusters_of_weight[1:min(2, length(clusters_of_weight))]  # Max 2 per weight
            if example_count >= 5
                break
            end
            
            example_count += 1
            
            total_loops = sum(values(cluster.multiplicities))
            loop_details = String[]
            for (loop_id, mult) in cluster.multiplicities
                push!(loop_details, "loop_$loop_id appears $mult times")
            end
            
            println("  Example $example_count: Weight $(cluster.weight)")
            println("    Composition: $(join(loop_details, "; "))")
            println("    Total loop instances: $total_loops")
        end
    end
end

function test_smaller_cases()
    """Test on smaller cases to verify correctness."""
    println("\n" * "="^50)
    println("Testing smaller cases for verification")
    println("="^50)
    
    test_cases = [
        (3, 1, 4, 4),  # 3x3 lattice, site 1, max weight 4
        (4, 1, 5, 4),  # 4x4 lattice, site 1, max weight 5
        (5, 1, 6, 5),  # 5x5 lattice, site 1, max weight 6
    ]
    
    results = Tuple{Int, Int}[]
    
    for (L, site, max_cluster_weight, max_loop_weight) in test_cases
        println("\nTesting $(L)x$(L) lattice, max weight $max_cluster_weight:")
        
        adj_matrix = create_periodic_square_lattice(L)
        enumerator = ClusterEnumerator(adj_matrix)
        
        start_time = time()
        try
            clusters = enumerate_connected_clusters(enumerator, site, max_cluster_weight, max_loop_weight)
            end_time = time()
            
            count = length(clusters)
            push!(results, (L, count))
            
            println("  Found $count clusters in $(round(end_time - start_time, digits=3))s")
            
        catch e
            end_time = time()
            println("  ERROR after $(round(end_time - start_time, digits=3))s: $e")
            push!(results, (L, 0))
        end
    end
    
    println("\nPattern analysis:")
    for (L, count) in results
        println("  $(L)x$(L): $count clusters")
    end
    
    return results
end

# function test_weight_8_vertex_degrees()
#     """Test if all weight-8 loops have all vertices with degree = 2."""
#     println("\n" * "="^60)
#     println("TESTING WEIGHT-8 LOOP VERTEX DEGREES")
#     println("="^60)
    
#     L = 10
#     adj_matrix = create_periodic_square_lattice(L)
#     enumerator = ClusterEnumerator(adj_matrix)
    
#     println("Parameters:")
#     println("  Lattice size: $(L)x$(L)")
#     println("  Testing weight-8 loops for degree constraint")
    
#     # Get all loops up to weight 8
#     println("\nEnumerating loops...")
#     start_time = time()
#     all_loops = find_all_loops_in_graph(enumerator, 8)
#     end_time = time()
    
#     # Filter weight-8 loops
#     weight_8_loops = [loop for loop in all_loops if loop.weight == 8]
#     println("Found $(length(weight_8_loops)) weight-8 loops in $(round(end_time - start_time, digits=3))s")
    
#     if isempty(weight_8_loops)
#         println("⚠️  No weight-8 loops found to test")
#         return true, 0, 0
#     end
    
#     all_degree_2 = true
#     degree_2_count = 0
#     higher_degree_count = 0
    
#     println("\nAnalyzing vertex degrees...")
#     for (i, loop) in enumerate(weight_8_loops)
#         # Calculate degree of each vertex in this loop
#         degree = Dict{Int, Int}()
#         for v in loop.vertices
#             degree[v] = 0
#         end
        
#         for (u, v) in loop.edges
#             degree[u] += 1
#             degree[v] += 1
#         end
        
#         # Check if all vertices have degree 2
#         loop_all_degree_2 = true
#         max_degree = 0
#         for v in loop.vertices
#             if degree[v] != 2
#                 loop_all_degree_2 = false
#                 all_degree_2 = false
#             end
#             max_degree = max(max_degree, degree[v])
#         end
        
#         if loop_all_degree_2
#             degree_2_count += 1
#         else
#             higher_degree_count += 1
#             if i <= 3  # Show first few examples
#                 println("  Loop $i: vertices=$(loop.vertices)")
#                 println("    Degrees: $([v => degree[v] for v in sort(loop.vertices)])")
#                 println("    Max degree: $max_degree")
#             end
#         end
#     end
    
#     println("\nRESULTS:")
#     println("Weight-8 loops with all degree = 2: $degree_2_count")
#     println("Weight-8 loops with some degree > 2: $higher_degree_count")
#     println("All weight-8 loops have degree = 2: $(all_degree_2 ? "✅ YES" : "❌ NO")")
    
#     if !all_degree_2
#         println("\n⚠️  ISSUE: Some weight-8 loops have vertices with degree > 2")
#         println("   This violates the loop definition (all vertices should have degree ≥ 2)")
#         println("   Loops are defined as connected subgraphs where every vertex has degree ≥ 2")
#         println("   For simple cycles, all vertices should have degree exactly 2")
#     end
    
#     return all_degree_2, degree_2_count, higher_degree_count
# end

function main()
    println("Julia Cluster Enumeration Test Suite")
    println("="^60)
    
    # Test smaller cases first
    test_smaller_cases()
    
    # # Test weight-8 vertex degrees
    # weight_8_passed, degree_2_count, higher_degree_count = test_weight_8_vertex_degrees()
    
    # Test the main 10x10 case
    println("\n" * "="^60)
    println("MAIN TEST: 10x10 lattice")
    println("="^60)
    
    total_clusters = test_10x10_weight_7_cluster_enumeration()
    
    println("\n" * "="^60)
    println("FINAL RESULTS")
    println("="^60)
    
    if total_clusters !== nothing
        println("10x10 lattice, max weight 7: $total_clusters connected clusters")
        println("Expected: 28")
        println("Weight 7 test: $(total_clusters == 28 ? "✅ PASS" : "❌ FAIL")")
    else
        println("Weight 7 test: ❌ FAIL (algorithm error)")
    end
    
    # println("Weight-8 degree test: $(weight_8_passed ? "✅ PASS" : "❌ FAIL")")
    # println("  Weight-8 loops with degree=2: $degree_2_count")
    # println("  Weight-8 loops with degree>2: $higher_degree_count")
    
    overall_success = (total_clusters == 28)
    println("\nOVERALL: $(overall_success ? "✅ ALL TESTS PASSED" : "❌ SOME TESTS FAILED")")
    
    return overall_success
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end