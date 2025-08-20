#!/usr/bin/env julia

"""
test_cluster_comprehensive.jl

Comprehensive test suite for cluster enumeration validation.
Combines all validation tests for the weight 10, size 11, PBC results.
"""

include("test_cluster_validation_simple.jl")
include("test_two_loop_clusters.jl") 
include("test_translation_symmetry_complete.jl")
include("test_translation_w6.jl")
include("test_deduplication.jl")

using Test

function run_all_cluster_tests()
    """Run the complete test suite."""
    println("🧪 Comprehensive Cluster Enumeration Test Suite")
    println("="^60)
    println("Testing cluster enumeration results with COMPLETE validations:")
    println("This validates ALL requirements for ALL clusters (no sampling):")
    println("  ✅ Multiplicity ≤ 2 (ALL clusters tested)")
    println("  ✅ Two-loop clusters have weights (4,4) or (4,6) (ALL tested)")
    println("  ✅ Two-loop clusters share 1, 2, or 4 vertices (ALL tested)")
    println("  ✅ Translation symmetry: 121 clusters per equivalence class (ALL tested)")
    println("  ✅ Translation symmetry on W6 data: Perfect equivalence classes")
    println("  ✅ Deduplication: Correct cluster counts (121 w4, 242 w6)")
    println("  ✅ Canonical forms: Proper geometric equivalence detection")
    println()
    
    # Run all test suites
    results = []
    
    println("🔍 Running Basic Validation Tests...")
    basic_result = run_basic_validation()
    push!(results, ("Basic Validation", basic_result))
    
    println("\n🔍 Running Two-Loop Cluster Analysis...")
    two_loop_result = analyze_two_loop_clusters()
    push!(results, ("Two-Loop Analysis", two_loop_result))
    
    println("\n🔍 Running Complete Translation Symmetry Test...")
    symmetry_result = test_translation_symmetry_complete()
    push!(results, ("Complete Translation Symmetry", symmetry_result))
    
    println("\n🔍 Running W6 Translation Symmetry Test...")
    w6_symmetry_result = test_translation_symmetry_w6_interface()
    push!(results, ("W6 Translation Symmetry", w6_symmetry_result))
    
    println("\n🔍 Running Deduplication Tests...")
    dedup_result = test_deduplication_functionality()
    push!(results, ("Deduplication", dedup_result))
    
    # Final summary
    println("\n" * "="^60)
    println("🏁 COMPREHENSIVE TEST SUITE SUMMARY")
    println("="^60)
    
    all_passed = true
    for (test_name, passed) in results
        status = passed ? "✅ PASS" : "❌ FAIL"
        println("$status - $test_name")
        all_passed = all_passed && passed
    end
    
    println("\n" * "="^60)
    if all_passed
        println("🎉 ALL TESTS PASSED!")
        println("✅ Cluster enumeration is fully validated")
        println("✅ All multiplicity constraints satisfied")
        println("✅ All loop weight constraints satisfied")
        println("✅ All vertex sharing constraints satisfied (ALL clusters)")
        println("✅ Translation symmetry preserved (ALL equivalence classes)")
        println("✅ Perfect translation symmetry on W6 data (121 per class)")
        println("✅ Deduplication working correctly")
        println("✅ Translation-aware canonical forms implemented")
        println("✅ COMPREHENSIVE testing completed - NO SAMPLING USED")
    else
        println("⚠️  SOME TESTS FAILED!")
        println("❌ Review the detailed results above")
    end
    
    return all_passed
end

# Run the comprehensive test if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Comprehensive Cluster Validation" begin
        @test run_all_cluster_tests()
    end
end