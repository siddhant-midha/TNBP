#!/usr/bin/env julia

"""
test_deduplication.jl

Test suite for cluster deduplication validation.
Verifies that the coordinate-based canonical form deduplication works correctly.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Serialization
using Test

include("test_utils.jl")

function load_test_enumeration_data()
    """Load the most recent test enumeration for validation."""
    data, filename = load_latest_test_file()
    return data
end

function test_cluster_count_expectations()
    """Test that cluster counts match theoretical expectations."""
    println("🔍 Testing Cluster Count Expectations")
    println("-"^40)
    
    data = load_test_enumeration_data()
    L = data.lattice_size
    max_weight = data.max_weight
    
    # Count clusters by weight
    weight_counts = Dict{Int, Int}()
    total_clusters = 0
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            weight_counts[cluster.weight] = get(weight_counts, cluster.weight, 0) + 1
            total_clusters += 1
        end
    end
    
    println("📊 Cluster distribution:")
    for weight in sort(collect(keys(weight_counts)))
        count = weight_counts[weight]
        println("  Weight $weight: $count clusters")
    end
    
    # Test expectations
    tests_passed = 0
    total_tests = 0
    
    # Test 1: Weight-4 clusters should be exactly 121 (one per site)
    total_tests += 1
    expected_w4 = L * L  # 121 for 11x11
    actual_w4 = get(weight_counts, 4, 0)
    if actual_w4 == expected_w4
        println("✅ Weight-4 clusters: $actual_w4 = $expected_w4 (expected)")
        tests_passed += 1
    else
        println("❌ Weight-4 clusters: $actual_w4 ≠ $expected_w4 (expected)")
    end
    
    # Test 2: Weight-6 clusters should be exactly 242 (2 per site)
    if max_weight >= 6
        total_tests += 1
        expected_w6 = 2 * L * L  # 242 for 11x11
        actual_w6 = get(weight_counts, 6, 0)
        if actual_w6 == expected_w6
            println("✅ Weight-6 clusters: $actual_w6 = $expected_w6 (expected)")
            tests_passed += 1
        else
            println("❌ Weight-6 clusters: $actual_w6 ≠ $expected_w6 (expected)")
        end
    end
    
    # Test 3: Each site should have the same number of clusters (translation symmetry)
    total_tests += 1
    clusters_per_site = [length(clusters) for clusters in values(data.clusters_by_site)]
    min_clusters = minimum(clusters_per_site)
    max_clusters = maximum(clusters_per_site)
    
    if min_clusters == max_clusters
        println("✅ Translation symmetry: All sites have $min_clusters clusters")
        tests_passed += 1
    else
        println("❌ Translation symmetry broken: sites have $min_clusters to $max_clusters clusters")
    end
    
    success_rate = tests_passed / total_tests
    println("📊 Test results: $tests_passed/$total_tests passed ($(round(success_rate*100))%)")
    
    return tests_passed == total_tests
end

function test_multiplicity_constraints()
    """Test that all multiplicities are ≤ 2 and structurally valid."""
    println("\n🔍 Testing Multiplicity Constraints")
    println("-"^40)
    
    data = load_test_enumeration_data()
    
    # Collect all clusters and their multiplicities
    all_clusters = []
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            push!(all_clusters, (site, cluster))
        end
    end
    
    # Check multiplicity constraints
    violations = []
    max_multiplicity = 0
    multiplicity_dist = Dict{Int, Int}()
    
    for (site, cluster) in all_clusters
        for (loop_id, multiplicity) in cluster.multiplicities
            max_multiplicity = max(max_multiplicity, multiplicity)
            multiplicity_dist[multiplicity] = get(multiplicity_dist, multiplicity, 0) + 1
            
            if multiplicity > 2
                push!(violations, (site, cluster.weight, multiplicity))
            end
        end
    end
    
    println("📊 Multiplicity distribution:")
    for mult in sort(collect(keys(multiplicity_dist)))
        count = multiplicity_dist[mult]
        println("  Multiplicity $mult: $count occurrences")
    end
    
    tests_passed = 0
    total_tests = 2
    
    # Test 1: No multiplicities > 2
    if isempty(violations)
        println("✅ All multiplicities ≤ 2 (max found: $max_multiplicity)")
        tests_passed += 1
    else
        println("❌ Found $(length(violations)) multiplicity violations")
        for (site, weight, mult) in violations[1:min(5, end)]
            println("    Site $site, weight $weight: multiplicity = $mult")
        end
    end
    
    # Test 2: For current test data, all multiplicities should be 1
    # (since we're testing weight ≤ 6 single-loop clusters)
    if max_multiplicity == 1
        println("✅ All multiplicities = 1 (expected for current test case)")
        tests_passed += 1
    else
        println("❌ Found multiplicities > 1 in single-loop test case")
    end
    
    success_rate = tests_passed / total_tests
    println("📊 Test results: $tests_passed/$total_tests passed ($(round(success_rate*100))%)")
    
    return tests_passed == total_tests
end

function test_canonical_form_consistency()
    """Test that canonical forms are working correctly for deduplication."""
    println("\n🔍 Testing Canonical Form Consistency")
    println("-"^40)
    
    data = load_test_enumeration_data()
    L = data.lattice_size
    
    # Get all weight-4 loops and test their canonical forms
    weight4_loops = [loop for loop in data.all_loops if loop.weight == 4]
    println("📊 Testing $(length(weight4_loops)) weight-4 loops")
    
    # Reproduce the canonical form function for testing
    function create_coordinate_canonical_form_test(loop, L::Int)
        function site_to_coords(site::Int)
            i = div(site - 1, L) + 1
            j = mod(site - 1, L) + 1
            return (i, j)
        end
        
        coords = [site_to_coords(v) for v in sort(loop.vertices)]
        
        # Check for periodic wrapping by examining coordinate spans
        i_values = [coord[1] for coord in coords]
        j_values = [coord[2] for coord in coords]
        
        i_span = maximum(i_values) - minimum(i_values)
        j_span = maximum(j_values) - minimum(j_values)
        
        # If span > L/2, the pattern likely wraps around periodic boundaries
        if i_span > L/2 || j_span > L/2
            unwrapped_coords = []
            
            for coord in coords
                i, j = coord
                if i_span > L/2
                    min_i = minimum(i_values)
                    if i - min_i > L/2
                        i = i - L
                    end
                end
                
                if j_span > L/2
                    min_j = minimum(j_values)
                    if j - min_j > L/2
                        j = j - L
                    end
                end
                
                push!(unwrapped_coords, (i, j))
            end
            
            coords = unwrapped_coords
        end
        
        # Normalize coordinates to start from (0, 0)
        min_i = minimum(coord[1] for coord in coords)
        min_j = minimum(coord[2] for coord in coords)
        
        normalized_coords = [(coord[1] - min_i, coord[2] - min_j) for coord in coords]
        normalized_coords = sort(normalized_coords)
        
        return (normalized_coords, loop.weight)
    end
    
    # Test canonical forms
    canonical_forms = Set()
    for loop in weight4_loops
        canonical = create_coordinate_canonical_form_test(loop, L)
        push!(canonical_forms, canonical)
    end
    
    tests_passed = 0
    total_tests = 2
    
    # Test 1: All weight-4 loops should have the same canonical form
    if length(canonical_forms) == 1
        println("✅ All weight-4 loops have the same canonical form")
        tests_passed += 1
    else
        println("❌ Found $(length(canonical_forms)) different canonical forms for weight-4 loops")
    end
    
    # Test 2: Similar test for weight-6 if available
    weight6_loops = [loop for loop in data.all_loops if loop.weight == 6]
    if !isempty(weight6_loops)
        println("📊 Testing $(length(weight6_loops)) weight-6 loops")
        
        w6_canonical_forms = Set()
        for loop in weight6_loops
            canonical = create_coordinate_canonical_form_test(loop, L)
            push!(w6_canonical_forms, canonical)
        end
        
        # Weight-6 loops should have exactly 2 canonical forms (2 different patterns)
        if length(w6_canonical_forms) == 2
            println("✅ Weight-6 loops have exactly 2 canonical forms")
            tests_passed += 1
        else
            println("❌ Found $(length(w6_canonical_forms)) canonical forms for weight-6 loops (expected 2)")
        end
    else
        println("ℹ️  No weight-6 loops to test")
        tests_passed += 1  # Skip this test
    end
    
    success_rate = tests_passed / total_tests
    println("📊 Test results: $tests_passed/$total_tests passed ($(round(success_rate*100))%)")
    
    return tests_passed == total_tests
end

function run_deduplication_tests()
    """Run all deduplication tests and return overall success."""
    println("🧪 Cluster Deduplication Test Suite")
    println("="^50)
    
    test_results = []
    
    # Run individual test suites
    push!(test_results, ("Cluster Count Expectations", test_cluster_count_expectations()))
    push!(test_results, ("Multiplicity Constraints", test_multiplicity_constraints()))
    push!(test_results, ("Canonical Form Consistency", test_canonical_form_consistency()))
    
    # Summary
    println("\n" * "="^50)
    println("🏁 DEDUPLICATION TEST SUMMARY")
    println("="^50)
    
    all_passed = true
    for (test_name, passed) in test_results
        status = passed ? "✅ PASS" : "❌ FAIL"
        println("$status - $test_name")
        all_passed = all_passed && passed
    end
    
    println("\n" * "="^50)
    if all_passed
        println("🎉 ALL DEDUPLICATION TESTS PASSED!")
        println("✅ Cluster counts match expectations")
        println("✅ Multiplicity constraints satisfied") 
        println("✅ Canonical forms working correctly")
        println("✅ Deduplication system is fully functional")
    else
        println("⚠️  SOME DEDUPLICATION TESTS FAILED!")
        println("❌ Review the detailed results above")
    end
    
    return all_passed
end

# Export for use in comprehensive test suite
function test_deduplication_functionality()
    """Public interface for the comprehensive test suite."""
    return run_deduplication_tests()
end

# Run tests if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Deduplication Tests" begin
        @test run_deduplication_tests()
    end
end