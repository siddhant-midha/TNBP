#!/usr/bin/env julia

"""
test_trace_redundancy.jl

Trace through the redundancy removal process to find the exact bug.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Serialization

function load_latest_cluster_file()
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    matching_files = filter(f -> contains(f, "L11") && contains(f, "w10") && contains(f, "periodic") && endswith(f, ".jld2"), files)
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end

function trace_redundancy_removal_bug()
    println("🔍 Tracing Redundancy Removal Bug")
    println("="^40)
    
    data = load_latest_cluster_file()
    
    # Get first few weight-4 clusters 
    weight4_clusters = []
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            if cluster.weight == 4
                push!(weight4_clusters, (site, cluster))
                if length(weight4_clusters) >= 10
                    break
                end
            end
        end
        if length(weight4_clusters) >= 10
            break
        end
    end
    
    println("📊 Analyzing first $(length(weight4_clusters)) weight-4 clusters:")
    
    # Check what canonical_cluster_signature returns for each
    signatures = []
    for (i, (site, cluster)) in enumerate(weight4_clusters)
        sig = canonical_cluster_signature(cluster)
        push!(signatures, sig)
        
        println("  $i. Site $site:")
        println("     Loop IDs: $(cluster.loop_ids)")
        println("     Multiplicities: $(cluster.multiplicities)")
        println("     canonical_cluster_signature: $sig")
    end
    
    # Check if signatures are actually unique
    unique_sigs = Set(signatures)
    println("\n📊 Signature analysis:")
    println("  Total signatures: $(length(signatures))")
    println("  Unique signatures: $(length(unique_sigs))")
    
    if length(unique_sigs) == length(signatures)
        println("❌ BUG FOUND: All signatures are unique - no deduplication possible!")
        println("   The canonical_cluster_signature function is NOT creating canonical forms")
    else
        println("✅ Some signatures are duplicated - deduplication should work")
    end
    
    # Show the unique signatures
    println("\n🔍 Unique signatures found:")
    for (i, sig) in enumerate(collect(unique_sigs))
        println("  $i: $sig")
    end
    
    # Test the remove_cluster_redundancies function directly
    println("\n🧪 Testing remove_cluster_redundancies directly...")
    test_clusters = [cluster for (site, cluster) in weight4_clusters]
    
    println("Input: $(length(test_clusters)) clusters")
    unique_result = remove_cluster_redundancies(test_clusters)
    println("Output: $(length(unique_result)) clusters")
    
    if length(unique_result) == length(test_clusters)
        println("❌ No deduplication occurred!")
    else
        println("✅ Deduplication worked: $(length(test_clusters) - length(unique_result)) clusters removed")
    end
    
    # The issue is clear: canonical_cluster_signature uses loop_id, which varies
    # We need to use canonical loop representations instead
    println("\n💡 SOLUTION NEEDED:")
    println("   canonical_cluster_signature should use canonical loop representations,")
    println("   not raw loop_ids, to create truly canonical cluster signatures.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    trace_redundancy_removal_bug()
end