#!/usr/bin/env julia

"""
test_final_dedup.jl

Test the final global deduplication logic to verify it works correctly.
This loads the 6x6 test results and analyzes the deduplication behavior.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("../functions/loop_enumeration_square.jl")

using Serialization

function load_test_file()
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    matching_files = filter(f -> contains(f, "test_global_dedup") && contains(f, "L6") && contains(f, "w4") && contains(f, "periodic") && endswith(f, ".jld2"), files)
    if isempty(matching_files)
        error("No test file found")
    end
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end

function test_canonical_signatures()
    println("🔍 Testing Canonical Signature Creation")
    println("="^50)
    
    data = load_test_file()
    
    # Get all weight-4 clusters
    weight4_clusters = []
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            if cluster.weight == 4
                push!(weight4_clusters, (site, cluster))
            end
        end
    end
    
    println("📊 Found $(length(weight4_clusters)) weight-4 clusters total")
    
    # Test both signature functions
    println("\n🧪 Testing both signature functions:")
    
    # Test the original canonical_cluster_signature
    original_sigs = Set()
    for (site, cluster) in weight4_clusters
        sig = canonical_cluster_signature(cluster)
        push!(original_sigs, sig)
    end
    
    println("  Original canonical_cluster_signature: $(length(original_sigs)) unique signatures")
    
    # Test the new create_canonical_cluster_signature  
    new_sigs = Set()
    for (site, cluster) in weight4_clusters
        sig = create_canonical_cluster_signature(cluster, data.all_loops)
        push!(new_sigs, sig)
    end
    
    println("  New create_canonical_cluster_signature: $(length(new_sigs)) unique signatures")
    
    # Show the difference
    if length(new_sigs) < length(original_sigs)
        println("  ✅ The new method finds fewer unique signatures - deduplication should work!")
        println("  Expected reduction: $(length(original_sigs)) → $(length(new_sigs))")
        
        # Show some example loops to understand the pattern
        println("\n📋 Example loop analysis (first 5 weight-4 clusters):")
        for (i, (site, cluster)) in enumerate(weight4_clusters[1:min(5, end)])
            println("  Cluster $i (site $site):")
            println("    Loop ID: $(cluster.loop_ids[1])")
            loop = data.all_loops[cluster.loop_ids[1]]
            canonical_loop = canonical_loop_representation(loop)
            println("    Loop vertices: $(loop.vertices)")
            println("    Canonical repr: $canonical_loop")
            
            original_sig = canonical_cluster_signature(cluster)
            new_sig = create_canonical_cluster_signature(cluster, data.all_loops)
            println("    Original sig: $original_sig")
            println("    New sig: $new_sig")
            println()
        end
        
    else
        println("  ❌ No improvement - both methods give same number of signatures")
    end
    
    # Test expected translation symmetry for 6x6 lattice
    if length(new_sigs) > 0
        expected_clusters_per_canonical = 36 ÷ length(new_sigs)
        total_expected = length(new_sigs) * expected_clusters_per_canonical
        println("\n📐 Translation symmetry analysis:")
        println("  Unique canonical forms: $(length(new_sigs))")
        println("  Expected clusters per form: $expected_clusters_per_canonical")
        println("  Expected total: $total_expected")
        println("  Actual total: $(length(weight4_clusters))")
        
        if total_expected == length(weight4_clusters)
            println("  ✅ Perfect translation symmetry!")
        else
            println("  ⚠️  Translation symmetry mismatch")
        end
    end
end

function test_manual_deduplication()
    println("\n🔧 Testing Manual Deduplication Process")
    println("="^50)
    
    data = load_test_file()
    
    # Manually apply the deduplication logic
    clusters_by_site = data.clusters_by_site
    all_loops = data.all_loops
    L = data.lattice_size
    
    # Collect all clusters with their sites (weight 4 only)
    weight4_clusters_with_sites = []
    for (site, clusters) in clusters_by_site
        for cluster in clusters
            if cluster.weight == 4
                push!(weight4_clusters_with_sites, (site, cluster))
            end
        end
    end
    
    println("  Original weight-4 clusters: $(length(weight4_clusters_with_sites))")
    
    # Group by canonical cluster signature (using loop canonical forms)
    canonical_groups = Dict()
    
    for (site, cluster) in weight4_clusters_with_sites
        canonical_sig = create_canonical_cluster_signature(cluster, all_loops)
        
        if !haskey(canonical_groups, canonical_sig)
            canonical_groups[canonical_sig] = []
        end
        push!(canonical_groups[canonical_sig], (site, cluster))
    end
    
    println("  Unique canonical forms: $(length(canonical_groups))")
    
    # Show group sizes
    group_sizes = [length(group) for group in values(canonical_groups)]
    println("  Group sizes: $(sort(group_sizes))")
    
    if length(canonical_groups) == 1
        println("  ✅ All weight-4 clusters have the same canonical form!")
        println("  This means the deduplication should reduce 144 → 36")
    else
        println("  Multiple canonical forms found")
        for (i, (canonical_sig, group)) in enumerate(canonical_groups)
            println("    Group $i: $(length(group)) clusters")
            println("      Canonical sig: $canonical_sig")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_canonical_signatures()
    test_manual_deduplication()
end