#!/usr/bin/env julia

"""
test_two_loop_clusters.jl

Focused test for finding and validating two-loop clusters with weights (4,4) and (4,6).
"""

# Include necessary modules
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("test_utils.jl")

using Test
using Serialization

function count_shared_vertices(loop1::Loop, loop2::Loop)
    """Count the number of shared vertices between two loops."""
    vertices1 = Set(loop1.vertices)
    vertices2 = Set(loop2.vertices)
    return length(intersect(vertices1, vertices2))
end

function analyze_two_loop_clusters()
    """Find and analyze clusters with two different loops."""
    println("🔍 Analyzing Two-Loop Clusters")
    println("="^40)
    
    # Load data
    local data
    try
        data, filename = load_latest_L11_w10_file()
        println("✅ Successfully loaded data")
        
        # Print summary
        total_clusters = sum(length(clusters) for clusters in values(data.clusters_by_site))
        println("📊 Total clusters: $total_clusters, Total loops: $(length(data.all_loops))")
        
    catch e
        println("❌ Failed to load cluster data: $e")
        return false
    end
    
    # Find clusters with two different loops
    two_loop_clusters = []
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            if length(cluster.loop_ids) == 2
                push!(two_loop_clusters, (site, cluster))
            end
        end
    end
    
    println("📊 Found $(length(two_loop_clusters)) clusters with two different loops")
    
    if isempty(two_loop_clusters)
        println("⚠️  No two-loop clusters found!")
        return false
    end
    
    # Analyze weight patterns
    weight_patterns = Dict{Vector{Int}, Int}()
    vertex_sharing_counts = Dict{Int, Int}()
    
    valid_count = 0
    invalid_weight_count = 0
    invalid_vertex_count = 0
    
    println("\n🔍 Analyzing weight patterns and vertex sharing for ALL two-loop clusters...")
    
    # Full analysis of ALL clusters with detailed output for first 20
    println("\\n📋 Detailed analysis of first 20 clusters:")
    
    for (i, (site, cluster)) in enumerate(two_loop_clusters)
        # Get loop weights
        weights = [data.all_loops[id].weight for id in cluster.loop_ids]
        sort!(weights)
        
        # Count weight patterns
        weight_patterns[weights] = get(weight_patterns, weights, 0) + 1
        
        # Get loops and count shared vertices
        loop1 = data.all_loops[cluster.loop_ids[1]]
        loop2 = data.all_loops[cluster.loop_ids[2]]
        shared = count_shared_vertices(loop1, loop2)
        
        # Count vertex sharing patterns
        vertex_sharing_counts[shared] = get(vertex_sharing_counts, shared, 0) + 1
        
        # Validation for ALL clusters
        weight_valid = (weights == [4, 4] || weights == [4, 6])
        vertex_valid = (shared in [1, 2, 4])
        
        if weight_valid && vertex_valid
            valid_count += 1
        end
        
        if !weight_valid
            invalid_weight_count += 1
        end
        
        if !vertex_valid
            invalid_vertex_count += 1
        end
        
        # Detailed output for first 20 only
        if i <= 20
            status_weight = weight_valid ? "✅" : "❌"
            status_vertex = vertex_valid ? "✅" : "❌"
            
            println("  $i. Site $site: weights $weights $status_weight, shared vertices: $shared $status_vertex")
            
            if i <= 3
                println("     Loop 1 vertices: $(loop1.vertices)")
                println("     Loop 2 vertices: $(loop2.vertices)")
                println("     Intersection: $(intersect(Set(loop1.vertices), Set(loop2.vertices)))")
            end
        end
    end
    
    if length(two_loop_clusters) > 20
        println("\\n📊 ... analyzed all $(length(two_loop_clusters)) two-loop clusters (showed details for first 20)")
    end
    
    # Summary statistics
    println("\\n📊 WEIGHT PATTERN SUMMARY:")
    for (pattern, count) in sort(collect(weight_patterns))
        valid = (pattern == [4, 4] || pattern == [4, 6])
        status = valid ? "✅" : "❌"
        println("  $pattern: $count clusters $status")
    end
    
    println("\\n📊 VERTEX SHARING SUMMARY:")
    for (shared, count) in sort(collect(vertex_sharing_counts))
        valid = (shared in [1, 2, 4])
        status = valid ? "✅" : "❌"
        println("  $shared shared vertices: $count clusters $status")
    end
    
    # Final validation
    expected_weights = Set([[4, 4], [4, 6]])
    found_weights = Set(collect(keys(weight_patterns)))
    
    expected_vertex_counts = Set([1, 2, 4])
    found_vertex_counts = Set(collect(keys(vertex_sharing_counts)))
    
    weight_test_passed = issubset(found_weights, expected_weights)
    vertex_test_passed = issubset(found_vertex_counts, expected_vertex_counts)
    
    println("\\n" * "="^50)
    println("🏁 TWO-LOOP CLUSTER VALIDATION")
    println("="^50)
    
    println("✅ Weight patterns: $(weight_test_passed ? "PASS" : "FAIL")")
    if !weight_test_passed
        unexpected = setdiff(found_weights, expected_weights)
        println("   Unexpected weights: $unexpected")
    end
    
    println("✅ Vertex sharing: $(vertex_test_passed ? "PASS" : "FAIL")")
    if !vertex_test_passed
        unexpected = setdiff(found_vertex_counts, expected_vertex_counts)
        println("   Unexpected vertex counts: $unexpected")
    end
    
    # Check for expected patterns
    has_4_4 = [4, 4] in keys(weight_patterns)
    has_4_6 = [4, 6] in keys(weight_patterns)
    
    println("✅ Found (4,4) clusters: $(has_4_4 ? "YES ($(weight_patterns[[4,4]]) clusters)" : "NO")")
    println("✅ Found (4,6) clusters: $(has_4_6 ? "YES ($(weight_patterns[[4,6]]) clusters)" : "NO")")
    
    overall_pass = weight_test_passed && vertex_test_passed && has_4_4 && has_4_6
    
    if overall_pass
        println("🎉 All two-loop cluster tests PASSED!")
    else
        println("⚠️  Some two-loop cluster tests FAILED!")
    end
    
    return overall_pass
end

# Run the test if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Two-Loop Cluster Analysis" begin
        @test analyze_two_loop_clusters()
    end
end