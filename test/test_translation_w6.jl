#!/usr/bin/env julia

"""
test_translation_w6.jl

Test translation symmetry specifically on L11_W6 data.
Focus on equivalence class sizes rather than total number.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("test_utils.jl")

using Test
using Serialization

function test_translation_symmetry_w6()
    """Test translation symmetry on L11_W6 data focusing on equivalence class sizes."""
    println("🔍 Testing Translation Symmetry on L11_W6 Data")
    println("="^50)
    
    # Load L11 W6 data specifically
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    
    # Look for L11 W6 files, preferring the fixed ones
    w6_files = filter(f -> contains(f, "L11") && contains(f, "w6") && contains(f, "periodic") && endswith(f, ".jld2"), files)
    
    if isempty(w6_files)
        println("❌ No L11 W6 files found!")
        return false
    end
    
    # Prioritize fixed files
    fixed_w6_files = filter(f -> contains(f, "fixed"), w6_files)
    if !isempty(fixed_w6_files)
        latest_file = sort(fixed_w6_files)[end]
        println("🔧 Using fixed W6 data")
    else
        latest_file = sort(w6_files)[end]
        println("⚠️  Using unfixed W6 data")
    end
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading W6 data: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    data = loaded_data["data"]
    
    # Print summary
    total_clusters = sum(length(clusters) for clusters in values(data.clusters_by_site))
    L = data.lattice_size
    expected_class_size = L * L  # 11² = 121
    
    println("✅ Successfully loaded data")
    println("📊 Total clusters: $total_clusters, Lattice: $(L)×$(L)")
    println("📊 Expected equivalence class size: $expected_class_size")
    
    # Test 1: Uniform distribution
    println("\n🔍 Test 1: Checking uniform distribution across all sites...")
    
    site_cluster_counts = []
    for site in 1:(L*L)
        site_clusters = get(data.clusters_by_site, site, [])
        count = length(site_clusters)
        push!(site_cluster_counts, count)
    end
    
    min_count = minimum(site_cluster_counts)
    max_count = maximum(site_cluster_counts)
    avg_count = sum(site_cluster_counts) / length(site_cluster_counts)
    
    println("📊 Cluster counts per site:")
    println("  Min: $min_count, Max: $max_count, Average: $(round(avg_count, digits=2))")
    
    uniform_distribution = (min_count == max_count)
    println("✅ Perfect uniform distribution: $(uniform_distribution ? "YES" : "NO")")
    
    # Test 2: Equivalence class sizes - THE KEY TEST
    println("\n🔍 Test 2: Testing equivalence class sizes (THE KEY TEST)...")
    
    # Group clusters by canonical signature
    signature_to_clusters = Dict()
    
    println("📊 Computing canonical signatures for all clusters...")
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            # Use translation-aware signature instead of basic canonical signature
            signature = translation_aware_cluster_signature(cluster, data.all_loops, L)
            
            if !haskey(signature_to_clusters, signature)
                signature_to_clusters[signature] = []
            end
            push!(signature_to_clusters[signature], (site, cluster))
        end
    end
    
    n_equivalence_classes = length(signature_to_clusters)
    println("📊 Found $n_equivalence_classes unique equivalence classes")
    
    # Check equivalence class sizes - THIS IS THE CRUCIAL TEST
    perfect_classes = 0
    imperfect_classes = 0
    class_sizes = []
    size_violations = []
    
    for (signature, cluster_list) in signature_to_clusters
        class_size = length(cluster_list)
        push!(class_sizes, class_size)
        
        if class_size == expected_class_size
            perfect_classes += 1
        else
            imperfect_classes += 1
            if length(size_violations) < 10  # Show first 10 violations
                push!(size_violations, class_size)
            end
        end
    end
    
    println("📊 EQUIVALENCE CLASS SIZE ANALYSIS:")
    println("  Perfect classes (size $expected_class_size): $perfect_classes")
    println("  Imperfect classes: $imperfect_classes")
    
    if imperfect_classes > 0
        println("  Class size distribution: min=$(minimum(class_sizes)), max=$(maximum(class_sizes))")
        if !isempty(size_violations)
            println("  First 10 violation sizes: $size_violations")
        end
    end
    
    # Calculate success rate
    success_rate = perfect_classes / n_equivalence_classes
    println("  Success rate: $(round(success_rate * 100, digits=1))% of classes have correct size")
    
    # Test 3: Weight-specific analysis
    println("\n🔍 Test 3: Weight-specific equivalence class analysis...")
    
    # Group by weight first, then by canonical signature
    weight_class_analysis = Dict()
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            weight = cluster.weight
            # Use translation-aware signature
            signature = translation_aware_cluster_signature(cluster, data.all_loops, L)
            
            if !haskey(weight_class_analysis, weight)
                weight_class_analysis[weight] = Dict()
            end
            
            if !haskey(weight_class_analysis[weight], signature)
                weight_class_analysis[weight][signature] = 0
            end
            weight_class_analysis[weight][signature] += 1
        end
    end
    
    for weight in sort(collect(keys(weight_class_analysis)))
        weight_signatures = weight_class_analysis[weight]
        n_classes = length(weight_signatures)
        perfect_w = 0
        imperfect_w = 0
        
        for (signature, count) in weight_signatures
            if count == expected_class_size
                perfect_w += 1
            else
                imperfect_w += 1
            end
        end
        
        success_rate_w = perfect_w / n_classes
        println("📊 Weight $weight: $n_classes classes, $(perfect_w) perfect, $(imperfect_w) imperfect ($(round(success_rate_w * 100, digits=1))% success)")
    end
    
    # Final verdict
    all_classes_perfect = (imperfect_classes == 0)
    
    println("\n" * "="^50)
    println("🏁 TRANSLATION SYMMETRY TEST SUMMARY (W6)")
    println("="^50)
    
    println("✅ Uniform distribution: $(uniform_distribution ? "PASS" : "FAIL")")
    println("✅ Perfect equivalence class sizes: $(all_classes_perfect ? "PASS" : "FAIL")")
    println("📊 Overall success rate: $(round(success_rate * 100, digits=1))%")
    
    if all_classes_perfect
        println("🎉 PERFECT! All $n_equivalence_classes equivalence classes have exactly $expected_class_size members")
        println("✅ Translation symmetry is completely satisfied")
    else
        println("⚠️  Some equivalence classes don't have the expected size")
        println("📊 $perfect_classes/$n_equivalence_classes classes are perfect")
    end
    
    return all_classes_perfect
end

# Export function for comprehensive test suite
function test_translation_symmetry_w6_interface()
    """Public interface for comprehensive test suite."""
    return test_translation_symmetry_w6()
end

if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Translation Symmetry W6 Test" begin
        @test test_translation_symmetry_w6()
    end
end