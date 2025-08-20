#!/usr/bin/env julia

"""
investigate_cluster_structure.jl

Investigate the actual structure of weight-4 clusters to understand 
what's causing the 484 vs 121 discrepancy.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Serialization

function load_latest_11x11_file()
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    matching_files = filter(f -> contains(f, "L11") && contains(f, "w4") && contains(f, "periodic") && endswith(f, ".jld2"), files)
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end

function analyze_weight4_cluster_structure()
    println("🔍 Investigating Weight-4 Cluster Structure")
    println("="^50)
    
    data = load_latest_11x11_file()
    
    # Collect all weight-4 clusters
    weight4_clusters = []
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            if cluster.weight == 4
                push!(weight4_clusters, (site, cluster))
            end
        end
    end
    
    println("📊 Found $(length(weight4_clusters)) weight-4 clusters total")
    
    # Analyze cluster structures
    structure_counts = Dict()
    
    for (site, cluster) in weight4_clusters
        # Analyze the structure
        n_loops = length(cluster.loop_ids)
        multiplicities = sort(collect(values(cluster.multiplicities)))
        
        # Get loop weights
        loop_weights = []
        for loop_id in cluster.loop_ids
            loop = data.all_loops[loop_id]
            push!(loop_weights, loop.weight)
        end
        sort!(loop_weights)
        
        structure_key = (n_loops, multiplicities, loop_weights)
        structure_counts[structure_key] = get(structure_counts, structure_key, 0) + 1
    end
    
    println("\n📋 Weight-4 cluster structures:")
    for (structure, count) in sort(collect(structure_counts), by=x->x[2], rev=true)
        n_loops, multiplicities, loop_weights = structure
        println("  $count clusters: $n_loops loops, multiplicities=$multiplicities, weights=$loop_weights")
        
        # Verify weight calculation
        total_weight = sum(loop_weights[i] * multiplicities[i] for i in 1:length(loop_weights))
        if total_weight != 4
            println("    ❌ ERROR: Total weight = $total_weight ≠ 4")
        end
    end
    
    # Look at a few examples in detail
    println("\n🔍 Detailed examples (first 5 clusters):")
    for (i, (site, cluster)) in enumerate(weight4_clusters[1:min(5, end)])
        println("  Example $i (site $site):")
        println("    Loop IDs: $(cluster.loop_ids)")
        println("    Multiplicities: $(cluster.multiplicities)")
        println("    Total loops in cluster: $(cluster.total_loops)")
        println("    Weight: $(cluster.weight)")
        
        # Show loop details
        for loop_id in cluster.loop_ids
            loop = data.all_loops[loop_id]
            mult = cluster.multiplicities[loop_id]
            println("      Loop $loop_id (×$mult): weight=$(loop.weight), vertices=$(loop.vertices)")
        end
        println()
    end
    
    # Check if there are actually different canonical forms
    println("🧪 Testing canonical signatures:")
    signatures = Set()
    canonical_sigs = Set()
    
    for (site, cluster) in weight4_clusters[1:min(20, end)]  # Test first 20
        # Original signature (uses loop IDs)
        orig_sig = canonical_cluster_signature(cluster)
        push!(signatures, orig_sig)
        
        # Try to create a canonical signature based on loop content
        canonical_loops = []
        for loop_id in cluster.loop_ids
            loop = data.all_loops[loop_id]
            canonical_loop = canonical_loop_representation(loop)
            mult = cluster.multiplicities[loop_id]
            push!(canonical_loops, (canonical_loop, mult))
        end
        sort!(canonical_loops)
        canonical_sig = (tuple(canonical_loops...), cluster.weight)
        push!(canonical_sigs, canonical_sig)
    end
    
    println("  Original signatures (first 20): $(length(signatures)) unique")
    println("  Canonical signatures (first 20): $(length(canonical_sigs)) unique")
    
    if length(canonical_sigs) < length(signatures)
        println("  ✅ Canonical approach should work for deduplication!")
    else
        println("  ❌ No improvement with canonical approach")
    end
    
    # Test with all clusters if sample shows promise
    if length(canonical_sigs) < length(signatures)
        println("\n🔄 Testing all $(length(weight4_clusters)) clusters...")
        all_canonical_sigs = Set()
        
        for (site, cluster) in weight4_clusters
            canonical_loops = []
            for loop_id in cluster.loop_ids
                loop = data.all_loops[loop_id]
                canonical_loop = canonical_loop_representation(loop)
                mult = cluster.multiplicities[loop_id]
                push!(canonical_loops, (canonical_loop, mult))
            end
            sort!(canonical_loops)
            canonical_sig = (tuple(canonical_loops...), cluster.weight)
            push!(all_canonical_sigs, canonical_sig)
        end
        
        println("  All canonical signatures: $(length(all_canonical_sigs)) unique")
        println("  Expected reduction: 484 → $(length(all_canonical_sigs))")
        
        if length(all_canonical_sigs) == 121
            println("  🎯 Perfect! This matches the expected 121 unique forms")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_weight4_cluster_structure()
end