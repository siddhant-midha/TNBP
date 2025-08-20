#!/usr/bin/env julia

"""
test_redundancy_bug.jl

Investigate why redundancy removal is failing for weight-4 clusters.
There should be 121 weight-4 clusters, not 484.
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

function investigate_weight4_redundancy_bug()
    println("🐛 Investigating Weight-4 Redundancy Bug")
    println("="^50)
    
    data = load_latest_cluster_file()
    L = data.lattice_size
    
    # Find all weight-4 clusters (should be 121, not 484)
    weight4_clusters = []
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            if cluster.weight == 4
                push!(weight4_clusters, (site, cluster))
            end
        end
    end
    
    println("📊 Found $(length(weight4_clusters)) weight-4 clusters")
    println("📊 Expected: 121 (one per translation)")
    println("📊 Problem ratio: $(length(weight4_clusters) / 121)")
    
    # Analyze the structure of these clusters
    println("\n🔍 Analyzing cluster structures...")
    
    single_loop_count = 0
    multi_loop_count = 0
    mult2_count = 0
    
    cluster_signatures = Set{Tuple}()
    canonical_signatures = Set{Tuple}()
    
    for (site, cluster) in weight4_clusters[1:min(20, length(weight4_clusters))]  # Sample first 20
        signature = canonical_cluster_signature(cluster)
        push!(cluster_signatures, signature)
        
        # Try to find canonical form under translation
        # This is where the bug likely is
        canonical_cluster = canonical_cluster_under_translation(cluster, L, data.all_loops)
        canonical_sig = canonical_cluster_signature(canonical_cluster)
        push!(canonical_signatures, canonical_sig)
        
        println("  Site $site:")
        println("    Loop IDs: $(cluster.loop_ids)")
        println("    Multiplicities: $(cluster.multiplicities)")
        println("    Signature: $signature")
        println("    Canonical signature: $canonical_sig")
        
        if length(cluster.loop_ids) == 1
            single_loop_count += 1
            loop_id = cluster.loop_ids[1]
            if cluster.multiplicities[loop_id] == 2
                mult2_count += 1
                println("    ⚠️  Single loop with multiplicity 2!")
            end
        else
            multi_loop_count += 1
        end
        println()
    end
    
    println("📊 Sample analysis:")
    println("  Single-loop clusters: $single_loop_count")
    println("  Multi-loop clusters: $multi_loop_count") 
    println("  Multiplicity-2 clusters: $mult2_count")
    println("  Unique signatures in sample: $(length(cluster_signatures))")
    println("  Unique canonical signatures in sample: $(length(canonical_signatures))")
    
    # Check if canonical form finding is working
    if length(canonical_signatures) < length(cluster_signatures)
        println("✅ Canonical form reduction is working in sample")
    else
        println("❌ Canonical form reduction NOT working - each cluster has unique canonical form")
    end
    
    # Let's check a specific example: find clusters that should be equivalent
    println("\n🔍 Checking for equivalent clusters that should be merged...")
    
    # Take the first weight-4 cluster and see if we can find its translations
    if !isempty(weight4_clusters)
        ref_site, ref_cluster = weight4_clusters[1]
        ref_canonical = canonical_cluster_under_translation(ref_cluster, L, data.all_loops)
        ref_canonical_sig = canonical_cluster_signature(ref_canonical)
        
        println("Reference cluster (site $ref_site):")
        println("  Loop IDs: $(ref_cluster.loop_ids)")
        println("  Canonical signature: $ref_canonical_sig")
        
        # Count how many clusters have the same canonical signature
        matching_count = 0
        for (site, cluster) in weight4_clusters
            canonical_cluster = canonical_cluster_under_translation(cluster, L, data.all_loops)
            canonical_sig = canonical_cluster_signature(canonical_cluster)
            
            if canonical_sig == ref_canonical_sig
                matching_count += 1
                if matching_count <= 5  # Show first 5 matches
                    println("  Match $matching_count: Site $site, Loop IDs: $(cluster.loop_ids)")
                end
            end
        end
        
        println("📊 Found $matching_count clusters with same canonical form")
        println("📊 Expected: 121 (if this is the only canonical form)")
        
        if matching_count == 484
            println("❌ BUG CONFIRMED: All 484 clusters have the same canonical form!")
            println("   This means redundancy removal completely failed.")
        elseif matching_count == 121
            println("✅ This canonical form appears correctly 121 times")
            println("   The bug might be that there are multiple canonical forms when there should be one")
        else
            println("❓ Unexpected count - needs further investigation")
        end
    end
    
    return length(weight4_clusters)
end

# Translation utilities for debugging
function site_to_coords(site::Int, L::Int)
    i = div(site - 1, L) + 1
    j = mod(site - 1, L) + 1
    return (i, j)
end

function translate_site(site::Int, di::Int, dj::Int, L::Int)
    i, j = site_to_coords(site, L)
    new_i = mod(i - 1 + di, L) + 1
    new_j = mod(j - 1 + dj, L) + 1
    return (new_i - 1) * L + new_j
end

function translate_cluster(cluster::Cluster, di::Int, dj::Int, L::Int, all_loops::Vector{Loop})
    translated_loop_ids = Int[]
    translated_multiplicities = Dict{Int, Int}()
    
    for (loop_id, multiplicity) in cluster.multiplicities
        original_loop = all_loops[loop_id]
        translated_vertices = [translate_site(v, di, dj, L) for v in original_loop.vertices]
        translated_edges = [(translate_site(u, di, dj, L), translate_site(v, di, dj, L)) for (u, v) in original_loop.edges]
        translated_loop = Loop(sort(translated_vertices), sort([(min(u,v), max(u,v)) for (u,v) in translated_edges]), original_loop.weight)
        translated_canonical = canonical_loop_representation(translated_loop)
        
        translated_loop_id = -1
        for (i, loop) in enumerate(all_loops)
            if canonical_loop_representation(loop) == translated_canonical
                translated_loop_id = i
                break
            end
        end
        
        if translated_loop_id == -1
            error("Could not find translated loop")
        end
        
        push!(translated_loop_ids, translated_loop_id)
        translated_multiplicities[translated_loop_id] = multiplicity
    end
    
    return Cluster(sort(unique(translated_loop_ids)), translated_multiplicities, cluster.weight, cluster.total_loops)
end

function canonical_cluster_under_translation(cluster::Cluster, L::Int, all_loops::Vector{Loop})
    canonical_cluster = cluster
    canonical_signature = canonical_cluster_signature(cluster)
    
    for di in 0:(L-1)
        for dj in 0:(L-1)
            if di == 0 && dj == 0
                continue
            end
            
            translated_cluster = translate_cluster(cluster, di, dj, L, all_loops)
            translated_signature = canonical_cluster_signature(translated_cluster)
            
            if translated_signature < canonical_signature
                canonical_cluster = translated_cluster
                canonical_signature = translated_signature
            end
        end
    end
    
    return canonical_cluster
end

if abspath(PROGRAM_FILE) == @__FILE__
    investigate_weight4_redundancy_bug()
end