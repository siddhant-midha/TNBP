#!/usr/bin/env julia

"""
test_cluster_comprehensive.jl

Comprehensive test suite for cluster enumeration validation.
Combines all validation tests for the weight 10, size 11, PBC results.
"""

include("test_cluster_validation_simple.jl")
include("test_two_loop_clusters.jl") 
include("test_translation_symmetry.jl")

using Test

function run_all_cluster_tests()
    """Run the complete test suite."""
    println("🧪 Comprehensive Cluster Enumeration Test Suite")
    println("="^60)
    println("Testing weight 10, size 11, PBC cluster enumeration results")
    println("This validates all requirements:")
    println("  ✅ Multiplicity ≤ 2")
    println("  ✅ Two-loop clusters have weights (4,4) or (4,6)")
    println("  ✅ Two-loop clusters share 1, 2, or 4 vertices")
    println("  ✅ Translation symmetry: 121 clusters per equivalence class")
    println()
    
    # Run all test suites
    results = []
    
    println("🔍 Running Basic Validation Tests...")
    basic_result = run_basic_validation()
    push!(results, ("Basic Validation", basic_result))
    
    println("\n🔍 Running Two-Loop Cluster Analysis...")
    two_loop_result = analyze_two_loop_clusters()
    push!(results, ("Two-Loop Analysis", two_loop_result))
    
    println("\n🔍 Running Translation Symmetry Test...")
    symmetry_result = test_translation_symmetry_sample()
    push!(results, ("Translation Symmetry", symmetry_result))
    
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
        println("✅ All vertex sharing constraints satisfied")
        println("✅ Translation symmetry preserved")
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