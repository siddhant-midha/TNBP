#!/usr/bin/env julia

"""
test_weight4_cluster_sanity.jl

Sanity check: verify that weight-4 clusters are only single loops (4-cycles).
This test loads the saved cluster data and examines all weight-4 clusters.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
using Serialization

function test_weight4_clusters_are_single_loops()
    """
    Test that all weight-4 clusters consist of exactly one loop with weight 4.
    """
    
    println("🔍 Testing Weight-4 Cluster Sanity Check")
    println("="^60)
    
    # Load saved cluster data
    cluster_file = "../saved_clusters/periodic_clusters_L11_w10_2025-08-14T20-25-54-769.jld2"
    
    println("📂 Loading cluster data from: $(basename(cluster_file))")
    
    loaded_data = open(cluster_file, "r") do io
        deserialize(io)
    end
    
    data = loaded_data["data"]
    summary = loaded_data["summary"]
    
    println("✅ Loaded successfully!")
    println("   Lattice size: $(summary["lattice_size"])×$(summary["lattice_size"])")
    println("   Max weight: $(summary["max_weight"])")
    println("   Total loops: $(summary["total_loops"])")
    println("   Total clusters: $(summary["total_clusters"])")
    println()
    
    # Get clusters from saved data
    if hasfield(typeof(data), :clusters_by_site)
        # All-sites enumeration data - use site 1 clusters (PBC makes all equivalent)
        all_clusters = data.clusters_by_site[1]
    else
        # Single-site enumeration data
        all_clusters = data.clusters
    end
    
    all_loops = data.all_loops
    
    # Filter for weight-4 clusters
    weight4_clusters = [c for c in all_clusters if c.weight == 4]
    
    println("📊 Found $(length(weight4_clusters)) clusters with weight 4")
    println()
    
    # Analyze each weight-4 cluster
    all_valid = true
    
    for (i, cluster) in enumerate(weight4_clusters)
        println("Cluster $i:")
        println("  Total loops: $(cluster.total_loops)")
        println("  Loop IDs: $(cluster.loop_ids)")
        println("  Multiplicities: $(cluster.multiplicities)")
        println("  Weight: $(cluster.weight)")
        
        # Check if this cluster has exactly one loop
        if cluster.total_loops != 1
            println("  ❌ ERROR: Expected 1 loop, found $(cluster.total_loops)")
            all_valid = false
        end
        
        # Check if that single loop has weight 4
        if length(cluster.loop_ids) != 1
            println("  ❌ ERROR: Expected 1 loop ID, found $(length(cluster.loop_ids))")
            all_valid = false
        else
            loop_id = cluster.loop_ids[1]
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            
            println("  Loop details:")
            println("    Loop ID: $loop_id")
            println("    Multiplicity: $multiplicity")
            println("    Loop weight: $(loop.weight)")
            println("    Loop vertices: $(loop.vertices)")
            println("    Loop edges: $(loop.edges)")
            
            if multiplicity != 1
                println("  ❌ ERROR: Expected multiplicity 1, found $multiplicity")
                all_valid = false
            end
            
            if loop.weight != 4
                println("  ❌ ERROR: Expected loop weight 4, found $(loop.weight)")
                all_valid = false
            end
            
            if length(loop.edges) != 4
                println("  ❌ ERROR: Expected 4 edges, found $(length(loop.edges))")
                all_valid = false
            end
            
            if length(loop.vertices) != 4
                println("  ❌ ERROR: Expected 4 vertices, found $(length(loop.vertices))")
                all_valid = false
            end
            
            if all_valid
                println("  ✅ Valid: Single 4-cycle as expected")
            end
        end
        
        println()
        
        # Only show first few for brevity
        if i >= 5
            println("  ... (showing first 5 out of $(length(weight4_clusters)) total)")
            break
        end
    end
    
    # Summary
    println("="^60)
    if all_valid
        println("✅ PASS: All weight-4 clusters are single 4-cycles")
        println("   This confirms the cluster enumeration is working correctly.")
    else
        println("❌ FAIL: Found invalid weight-4 clusters")
        println("   The cluster enumeration may have issues.")
    end
    
    # Additional check: count unique 4-cycles
    unique_loops = Set()
    for cluster in weight4_clusters
        if cluster.total_loops == 1 && length(cluster.loop_ids) == 1
            loop_id = cluster.loop_ids[1]
            loop = all_loops[loop_id]
            if loop.weight == 4
                # Create canonical representation of the loop
                sorted_edges = sort([(min(e[1], e[2]), max(e[1], e[2])) for e in loop.edges])
                push!(unique_loops, sorted_edges)
            end
        end
    end
    
    println()
    println("📈 Statistics:")
    println("   Weight-4 clusters found: $(length(weight4_clusters))")
    println("   Unique 4-cycles represented: $(length(unique_loops))")
    
    # Check for expected number of 4-cycles on 11x11 lattice
    L = summary["lattice_size"]
    # Each site can be the top-left corner of a 2x2 square (4-cycle)
    expected_4cycles = L^2  # Each site has exactly one 4-cycle starting from it
    println("   Expected 4-cycles on $(L)×$(L) lattice: $expected_4cycles")
    
    if length(unique_loops) == expected_4cycles
        println("   ✅ Correct number of unique 4-cycles")
    else
        println("   ⚠️  Unexpected number of unique 4-cycles")
    end
    
    return all_valid
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    test_weight4_clusters_are_single_loops()
end