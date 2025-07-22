"""
test_cluster_enumeration.jl

Test suite for cluster enumeration in Julia.
"""

include("../functions/ClusterEnumeration.jl")

function test_10x10_cluster_enumeration()
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

function main()
    println("Julia Cluster Enumeration Test Suite")
    println("="^60)
    
    # Test smaller cases first
    test_smaller_cases()
    
    # Test the main 10x10 case
    println("\n" * "="^60)
    println("MAIN TEST: 10x10 lattice")
    println("="^60)
    
    total_clusters = test_10x10_cluster_enumeration()
    
    println("\n" * "="^60)
    println("FINAL RESULT")
    println("="^60)
    if total_clusters !== nothing
        println("10x10 lattice, max weight 7: $total_clusters connected clusters")
        println("Expected: 28")
        println("Match: $(total_clusters == 28 ? "✅ YES" : "❌ NO")")
        println("Algorithm completed successfully ✅")
        return total_clusters == 28
    else
        println("Algorithm failed ❌")
        return false
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end