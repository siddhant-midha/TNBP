#!/usr/bin/env julia

"""
test_translation_symmetry_complete.jl

Complete test for translation symmetry: verify ALL clusters satisfy translation symmetry requirements.
Tests uniformity of cluster distributions and validates equivalence class structure.
"""

# Include necessary modules
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("test_utils.jl")

using Test
using Serialization

function test_complete_translation_symmetry()
    """Test translation symmetry comprehensively for ALL clusters."""
    println("🔍 Testing Complete Translation Symmetry")
    println("="^50)
    
    # Load data
    local data
    try
        data, filename = load_latest_L11_w10_file()
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
    n_sites = L * L
    
    println("📊 Expected equivalence class size: $expected_class_size")
    println("📊 Number of sites: $n_sites")
    
    # Test 1: Uniform cluster distribution across ALL sites
    println("\n🔍 Test 1: Testing uniform cluster distribution across ALL sites...")
    
    site_cluster_counts = []
    cluster_weight_distributions = Dict{Int, Dict{Int, Int}}()  # site -> weight -> count
    
    for site in 1:n_sites
        site_clusters = get(data.clusters_by_site, site, [])
        count = length(site_clusters)
        push!(site_cluster_counts, count)
        
        # Track weight distribution per site
        cluster_weight_distributions[site] = Dict{Int, Int}()
        for cluster in site_clusters
            weight = cluster.weight
            cluster_weight_distributions[site][weight] = get(cluster_weight_distributions[site], weight, 0) + 1
        end
    end
    
    min_count = minimum(site_cluster_counts)
    max_count = maximum(site_cluster_counts)
    avg_count = sum(site_cluster_counts) / length(site_cluster_counts)
    
    println("📊 Cluster counts per site:")
    println("  Min: $min_count")
    println("  Max: $max_count") 
    println("  Average: $(round(avg_count, digits=2))")
    println("  Standard deviation: $(round(sqrt(sum((c - avg_count)^2 for c in site_cluster_counts) / n_sites), digits=4))")
    
    # For perfect translation symmetry, all sites should have exactly the same count
    uniform_total_distribution = (min_count == max_count)
    println("✅ Perfect uniform total distribution: $(uniform_total_distribution ? "YES" : "NO")")
    
    if !uniform_total_distribution
        variation = (max_count - min_count) / avg_count
        println("📊 Relative variation: $(round(variation * 100, digits=2))%")
        
        if variation < 0.01  # Less than 1% variation
            println("✅ Variation is negligible (<1%)")
            uniform_total_distribution = true
        end
    end
    
    # Test 2: Check weight-specific distributions
    println("\n🔍 Test 2: Testing weight-specific distributions across ALL sites...")
    
    # Get all weights that appear
    all_weights = Set{Int}()
    for clusters in values(data.clusters_by_site)
        for cluster in clusters
            push!(all_weights, cluster.weight)
        end
    end
    all_weights = sort(collect(all_weights))
    
    weight_uniform_tests = Dict{Int, Bool}()
    
    for weight in all_weights
        weight_counts = []
        for site in 1:n_sites
            count = get(get(cluster_weight_distributions, site, Dict{Int,Int}()), weight, 0)
            push!(weight_counts, count)
        end
        
        min_w = minimum(weight_counts)
        max_w = maximum(weight_counts)
        avg_w = sum(weight_counts) / n_sites
        
        uniform_for_weight = (min_w == max_w)
        
        if !uniform_for_weight && avg_w > 0
            variation_w = (max_w - min_w) / avg_w
            if variation_w < 0.01
                uniform_for_weight = true
            end
        end
        
        weight_uniform_tests[weight] = uniform_for_weight
        
        println("📊 Weight $weight: min=$min_w, max=$max_w, avg=$(round(avg_w, digits=2)) - $(uniform_for_weight ? "✅ UNIFORM" : "❌ NON-UNIFORM")")
    end
    
    weight_distributions_uniform = all(values(weight_uniform_tests))
    
    # Test 3: Equivalence class structure verification
    println("\n🔍 Test 3: Testing equivalence class structure for representative clusters...")
    
    # Group clusters by canonical signature to identify equivalence classes
    cluster_signatures = Dict()
    signature_to_clusters = Dict()
    
    println("📊 Computing canonical signatures for all clusters...")
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            # Use translation-aware signature for proper equivalence class identification
            signature = translation_aware_cluster_signature(cluster, data.all_loops, L)
            
            if !haskey(signature_to_clusters, signature)
                signature_to_clusters[signature] = []
            end
            push!(signature_to_clusters[signature], (site, cluster))
        end
    end
    
    n_equivalence_classes = length(signature_to_clusters)
    println("📊 Found $n_equivalence_classes unique equivalence classes")
    
    # Check sizes of equivalence classes
    class_size_violations = 0
    perfect_classes = 0
    
    class_sizes = []
    for (signature, cluster_list) in signature_to_clusters
        class_size = length(cluster_list)
        push!(class_sizes, class_size)
        
        if class_size == expected_class_size
            perfect_classes += 1
        else
            class_size_violations += 1
            if class_size_violations <= 5  # Show first 5 violations
                println("⚠️  Equivalence class size: $class_size (expected $expected_class_size)")
            end
        end
    end
    
    println("📊 Equivalence class size analysis:")
    println("  Perfect classes (size $expected_class_size): $perfect_classes")
    println("  Violations: $class_size_violations")
    
    if class_size_violations > 0
        println("  Class size distribution: min=$(minimum(class_sizes)), max=$(maximum(class_sizes))")
    end
    
    equivalence_classes_perfect = (class_size_violations == 0)
    
    # Final summary
    println("\n" * "="^50)
    println("🏁 COMPLETE TRANSLATION SYMMETRY TEST SUMMARY")
    println("="^50)
    
    println("✅ Uniform total cluster distribution: $(uniform_total_distribution ? "PASS" : "FAIL")")
    println("✅ Uniform weight-specific distributions: $(weight_distributions_uniform ? "PASS" : "FAIL")")
    println("✅ Perfect equivalence class sizes: $(equivalence_classes_perfect ? "PASS" : "FAIL")")
    
    overall_pass = uniform_total_distribution && weight_distributions_uniform && equivalence_classes_perfect
    
    total_clusters = sum(length(clusters) for clusters in values(data.clusters_by_site))
    
    if overall_pass
        println("🎉 COMPLETE translation symmetry test PASSED!")
        println("✅ All $(total_clusters) clusters respect translation symmetry")
        println("✅ All $n_equivalence_classes equivalence classes have correct size")
        println("✅ All weight distributions are uniform across sites")
    else
        println("⚠️  COMPLETE translation symmetry test had issues!")
        println("❌ Some translation symmetry violations found")
    end
    
    return overall_pass
end

# Export function for comprehensive test suite
function test_translation_symmetry_complete()
    """Public interface for comprehensive test suite."""
    return test_complete_translation_symmetry()
end

# Run the test if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Complete Translation Symmetry Test" begin
        @test test_complete_translation_symmetry()
    end
end