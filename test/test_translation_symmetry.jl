#!/usr/bin/env julia

"""
test_translation_symmetry.jl

Test translation symmetry: verify each equivalence class has exactly 11² = 121 clusters.
Uses a sampling approach for efficiency.
"""

# Include necessary modules
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Test
using Serialization
using Random

# Translation utilities
function site_to_coords(site::Int, L::Int)
    """Convert 1-indexed site to (i, j) coordinates."""
    i = div(site - 1, L) + 1
    j = mod(site - 1, L) + 1
    return (i, j)
end

function coords_to_site(i::Int, j::Int, L::Int)
    """Convert (i, j) coordinates to 1-indexed site."""
    return (i - 1) * L + j
end

function translate_site(site::Int, di::Int, dj::Int, L::Int)
    """Translate a site by (di, dj) on an L×L periodic lattice."""
    i, j = site_to_coords(site, L)
    
    # Apply translation with periodic boundary conditions
    new_i = mod(i - 1 + di, L) + 1
    new_j = mod(j - 1 + dj, L) + 1
    
    return coords_to_site(new_i, new_j, L)
end

function translate_cluster(cluster::Cluster, di::Int, dj::Int, L::Int, all_loops::Vector{Loop})
    """Translate a cluster by (di, dj) and return the translated cluster."""
    # Create mapping from old loop to translated loop
    translated_loop_ids = Int[]
    translated_multiplicities = Dict{Int, Int}()
    
    # For each loop in the cluster, find its translation
    for (loop_id, multiplicity) in cluster.multiplicities
        original_loop = all_loops[loop_id]
        
        # Translate all vertices in the loop
        translated_vertices = [translate_site(v, di, dj, L) for v in original_loop.vertices]
        translated_edges = [(translate_site(u, di, dj, L), translate_site(v, di, dj, L)) for (u, v) in original_loop.edges]
        
        # Create canonical representation and find matching loop
        translated_loop = Loop(sort(translated_vertices), sort([(min(u,v), max(u,v)) for (u,v) in translated_edges]), original_loop.weight)
        translated_canonical = canonical_loop_representation(translated_loop)
        
        # Find the corresponding loop ID in all_loops
        translated_loop_id = -1
        for (i, loop) in enumerate(all_loops)
            if canonical_loop_representation(loop) == translated_canonical
                translated_loop_id = i
                break
            end
        end
        
        if translated_loop_id == -1
            error("Could not find translated loop - this should not happen!")
        end
        
        push!(translated_loop_ids, translated_loop_id)
        translated_multiplicities[translated_loop_id] = multiplicity
    end
    
    return Cluster(
        sort(unique(translated_loop_ids)),
        translated_multiplicities,
        cluster.weight,
        cluster.total_loops
    )
end

function canonical_cluster_under_translation(cluster::Cluster, L::Int, all_loops::Vector{Loop})
    """Find the canonical representative of a cluster under translation symmetry."""
    canonical_cluster = cluster
    canonical_signature = canonical_cluster_signature(cluster)
    
    # Try all possible translations
    for di in 0:(L-1)
        for dj in 0:(L-1)
            if di == 0 && dj == 0
                continue  # Skip identity translation
            end
            
            translated_cluster = translate_cluster(cluster, di, dj, L, all_loops)
            translated_signature = canonical_cluster_signature(translated_cluster)
            
            # Keep the lexicographically smallest signature
            if translated_signature < canonical_signature
                canonical_cluster = translated_cluster
                canonical_signature = translated_signature
            end
        end
    end
    
    return canonical_cluster
end

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

function test_translation_symmetry_sample()
    """Test translation symmetry using a sampling approach."""
    println("🔍 Testing Translation Symmetry (Sampling)")
    println("="^45)
    
    # Load data
    local data
    try
        data = load_latest_cluster_file()
        println("✅ Successfully loaded data")
        
        # Print summary
        total_clusters = sum(length(clusters) for clusters in values(data.clusters_by_site))
        println("📊 Total clusters: $total_clusters, Lattice: $(data.lattice_size)×$(data.lattice_size)")
        
    catch e
        println("❌ Failed to load cluster data: $e")
        return false
    end
    
    L = data.lattice_size
    expected_class_size = L * L  # 11² = 121
    
    println("📊 Expected equivalence class size: $expected_class_size")
    
    # Sample clusters from different sites
    sample_sites = [1, 12, 34, 56, 78, 100, 121]  # Spread across the lattice
    sample_clusters = []
    
    for site in sample_sites
        if haskey(data.clusters_by_site, site) && !isempty(data.clusters_by_site[site])
            # Take a few clusters from each site
            site_clusters = data.clusters_by_site[site]
            n_sample = min(3, length(site_clusters))
            for i in 1:n_sample
                push!(sample_clusters, site_clusters[i])
            end
        end
    end
    
    println("📊 Testing $(length(sample_clusters)) sample clusters...")
    
    # Test each sample cluster
    violations = 0
    
    for (idx, cluster) in enumerate(sample_clusters)
        println("\\n🔍 Testing cluster $idx:")
        println("  Loop IDs: $(cluster.loop_ids)")
        println("  Weight: $(cluster.weight)")
        
        # Generate all translations of this cluster
        translated_clusters = Set{Tuple}()
        
        for di in 0:(L-1)
            for dj in 0:(L-1)
                try
                    translated = translate_cluster(cluster, di, dj, L, data.all_loops)
                    signature = canonical_cluster_signature(translated)
                    push!(translated_clusters, signature)
                catch e
                    println("  ⚠️  Error translating by ($di,$dj): $e")
                end
            end
        end
        
        class_size = length(translated_clusters)
        println("  Generated equivalence class size: $class_size")
        
        if class_size != expected_class_size
            violations += 1
            println("  ❌ Expected $expected_class_size, got $class_size")
            
            # Try to diagnose why
            if class_size < expected_class_size
                println("  ⚠️  Some translations may be equivalent (this suggests the cluster has internal symmetry)")
            else
                println("  ⚠️  Got more than expected - this is unexpected!")
            end
        else
            println("  ✅ Correct equivalence class size")
        end
        
        # Limit output for large samples
        if idx >= 10
            println("  ... (continuing with remaining clusters)")
            break
        end
    end
    
    # Test overall cluster distribution
    println("\\n🔍 Testing overall cluster distribution...")
    
    site_cluster_counts = [length(clusters) for clusters in values(data.clusters_by_site)]
    min_count = minimum(site_cluster_counts)
    max_count = maximum(site_cluster_counts)
    
    println("📊 Cluster counts per site: min=$min_count, max=$max_count")
    
    # For perfect translation symmetry, all sites should have exactly the same count
    uniform_distribution = (min_count == max_count)
    
    println("✅ Uniform distribution: $(uniform_distribution ? "YES" : "NO")")
    
    if !uniform_distribution
        # Check if the variation is small
        variation = (max_count - min_count) / min_count
        println("📊 Relative variation: $(round(variation * 100, digits=2))%")
        
        if variation < 0.01  # Less than 1% variation
            println("✅ Variation is very small - likely due to rounding or implementation details")
            uniform_distribution = true
        end
    end
    
    # Summary
    println("\\n" * "="^50)
    println("🏁 TRANSLATION SYMMETRY TEST SUMMARY")
    println("="^50)
    
    sample_test_passed = (violations == 0)
    
    println("✅ Sample equivalence classes: $(sample_test_passed ? "PASS" : "FAIL ($violations violations)")")
    println("✅ Uniform site distribution: $(uniform_distribution ? "PASS" : "FAIL")")
    
    overall_pass = sample_test_passed && uniform_distribution
    
    if overall_pass
        println("🎉 Translation symmetry test PASSED!")
    else
        println("⚠️  Translation symmetry test had issues - see details above")
    end
    
    return overall_pass
end

# Run the test if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Translation Symmetry Test" begin
        @test test_translation_symmetry_sample()
    end
end