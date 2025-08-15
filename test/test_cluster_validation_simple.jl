#!/usr/bin/env julia

"""
test_cluster_validation_simple.jl

Simplified validation test for cluster enumeration results.
Tests key properties without the full translation symmetry check.
"""

# Include necessary modules
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Test
using Serialization

function load_latest_cluster_file()
    """Load the most recent cluster enumeration file."""
    save_dir = "saved_clusters"
    
    if !isdir(save_dir)
        error("No saved_clusters directory found!")
    end
    
    # Look for files matching our criteria (weight 10, size 11, PBC)
    files = readdir(save_dir)
    
    # Filter for L11_w10 PBC files
    matching_files = filter(f -> contains(f, "L11") && contains(f, "w10") && contains(f, "periodic") && endswith(f, ".jld2"), files)
    
    if isempty(matching_files)
        error("No matching cluster files found for L=11, w=10, PBC!")
    end
    
    # Sort by timestamp and take the most recent
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    
    println("📖 Loading cluster data from: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    
    return loaded_data["data"]
end

function count_shared_vertices(loop1::Loop, loop2::Loop)
    """Count the number of shared vertices between two loops."""
    vertices1 = Set(loop1.vertices)
    vertices2 = Set(loop2.vertices)
    return length(intersect(vertices1, vertices2))
end

function run_basic_validation()
    """Run basic validation tests."""
    println("🧪 Basic Cluster Validation Test")
    println("="^50)
    
    # Load data
    local data
    try
        data = load_latest_cluster_file()
        println("✅ Successfully loaded data")
        
        # Print summary
        total_clusters = sum(length(clusters) for clusters in values(data.clusters_by_site))
        println("📊 Summary: $(data.lattice_size)×$(data.lattice_size) lattice, max weight $(data.max_weight)")
        println("📊 Total clusters: $total_clusters, Total loops: $(length(data.all_loops))")
        
    catch e
        println("❌ Failed to load cluster data: $e")
        return false
    end
    
    # Test 1: Multiplicity constraints
    println("\n🔍 Test 1: Checking multiplicity constraints...")
    mult_violations = 0
    mult2_count = 0
    sample_clusters = []
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            # Check multiplicity ≤ 2
            max_mult = maximum(values(cluster.multiplicities))
            if max_mult > 2
                mult_violations += 1
                if mult_violations <= 3
                    println("  ❌ Site $site: multiplicity $max_mult")
                end
            end
            
            # Count multiplicity=2 clusters for further analysis
            if any(mult == 2 for mult in values(cluster.multiplicities))
                mult2_count += 1
                if length(sample_clusters) < 10
                    push!(sample_clusters, (site, cluster))
                end
            end
        end
        
        if mult_violations > 10  # Early exit to avoid too much output
            break
        end
    end
    
    println("  📊 Multiplicity violations: $mult_violations")
    println("  📊 Clusters with multiplicity=2: $mult2_count")
    
    # Test 2: Sample analysis of multiplicity=2 clusters
    println("\n🔍 Test 2: Analyzing sample multiplicity=2 clusters...")
    weight_violations = 0
    vertex_violations = 0
    
    for (site, cluster) in sample_clusters[1:min(5, length(sample_clusters))]
        println("  📍 Site $site cluster:")
        println("    Loop IDs: $(cluster.loop_ids)")
        println("    Multiplicities: $(cluster.multiplicities)")
        
        # Check loop weights
        if length(cluster.loop_ids) == 1
            # Single loop with multiplicity 2
            loop_weight = data.all_loops[cluster.loop_ids[1]].weight
            println("    Single loop weight: $loop_weight")
            if loop_weight != 4
                weight_violations += 1
                println("    ❌ Expected weight 4, got $loop_weight")
            end
        elseif length(cluster.loop_ids) == 2
            # Two different loops
            weights = [data.all_loops[id].weight for id in cluster.loop_ids]
            sort!(weights)
            println("    Two loop weights: $weights")
            
            if !(weights == [4, 4] || weights == [4, 6])
                weight_violations += 1
                println("    ❌ Expected [4,4] or [4,6], got $weights")
            end
            
            # Check vertex sharing
            loop1 = data.all_loops[cluster.loop_ids[1]]
            loop2 = data.all_loops[cluster.loop_ids[2]]
            shared = count_shared_vertices(loop1, loop2)
            println("    Shared vertices: $shared")
            
            if !(shared in [1, 2, 4])
                vertex_violations += 1
                println("    ❌ Expected 1, 2, or 4 shared vertices, got $shared")
            end
        end
        println()
    end
    
    # Test 3: Translation symmetry sample check
    println("🔍 Test 3: Sample translation symmetry check...")
    
    # Pick a few sites and check if their cluster counts are similar
    # (not a full check, but gives us a sense)
    site_cluster_counts = [length(clusters) for clusters in values(data.clusters_by_site)]
    min_count = minimum(site_cluster_counts)
    max_count = maximum(site_cluster_counts)
    avg_count = sum(site_cluster_counts) / length(site_cluster_counts)
    
    println("  📊 Cluster counts per site: min=$min_count, max=$max_count, avg=$(round(avg_count, digits=1))")
    
    # For periodic boundary conditions, all sites should have similar cluster counts
    if max_count - min_count > avg_count * 0.1  # Allow 10% variation
        println("  ⚠️  Large variation in cluster counts - potential translation symmetry issue")
    else
        println("  ✅ Cluster counts are reasonably uniform across sites")
    end
    
    # Summary
    println("\n" * "="^50)
    println("🏁 VALIDATION SUMMARY")
    println("="^50)
    
    all_good = (mult_violations == 0 && weight_violations == 0 && vertex_violations == 0)
    
    println("✅ Multiplicity constraints: $(mult_violations == 0 ? "PASS" : "FAIL ($mult_violations violations)")")
    println("✅ Weight constraints: $(weight_violations == 0 ? "PASS" : "FAIL ($weight_violations violations)")")
    println("✅ Vertex sharing: $(vertex_violations == 0 ? "PASS" : "FAIL ($vertex_violations violations)")")
    
    if all_good
        println("🎉 Basic validation PASSED!")
    else
        println("⚠️  Some issues found - see details above")
    end
    
    return all_good
end

# Run the test if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Basic Cluster Validation" begin
        @test run_basic_validation()
    end
end