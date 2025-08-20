#!/usr/bin/env julia

"""
test_multiplicity_fix.jl

Test to verify that the multiplicity problem is fixed by checking recent enumeration results.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Serialization

function load_latest_test_file()
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    # Look for the most recent test files
    matching_files = filter(f -> contains(f, "test_") && contains(f, "L11") && endswith(f, ".jld2"), files)
    if isempty(matching_files)
        error("No test files found")
    end
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end

function test_multiplicity_constraints()
    println("🔍 Testing Multiplicity Constraints")
    println("="^50)
    
    data = load_latest_test_file()
    
    # Collect all clusters
    all_clusters = []
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            push!(all_clusters, (site, cluster))
        end
    end
    
    println("📊 Testing $(length(all_clusters)) total clusters")
    
    # Test multiplicity constraints
    multiplicity_violations = []
    max_multiplicity_found = 0
    multiplicity_distribution = Dict()
    
    for (site, cluster) in all_clusters
        for (loop_id, multiplicity) in cluster.multiplicities
            max_multiplicity_found = max(max_multiplicity_found, multiplicity)
            multiplicity_distribution[multiplicity] = get(multiplicity_distribution, multiplicity, 0) + 1
            
            # Check constraint: multiplicity should be ≤ 2
            if multiplicity > 2
                push!(multiplicity_violations, (site, cluster, loop_id, multiplicity))
            end
        end
    end
    
    println("📊 Multiplicity distribution:")
    for mult in sort(collect(keys(multiplicity_distribution)))
        count = multiplicity_distribution[mult]
        println("  Multiplicity $mult: $count occurrences")
    end
    
    println("📊 Maximum multiplicity found: $max_multiplicity_found")
    
    if isempty(multiplicity_violations)
        println("✅ All multiplicities ≤ 2: PASS")
        multiplicity_ok = true
    else
        println("❌ Found $(length(multiplicity_violations)) multiplicity violations:")
        for (site, cluster, loop_id, mult) in multiplicity_violations
            println("  Site $site, Loop $loop_id: multiplicity = $mult")
        end
        multiplicity_ok = false
    end
    
    # Test cluster structure constraints
    println("\n🔍 Testing cluster structure constraints...")
    
    structure_violations = []
    
    for (site, cluster) in all_clusters
        # Get loop weights and multiplicities
        loop_data = []
        for loop_id in cluster.loop_ids
            loop = data.all_loops[loop_id]
            mult = cluster.multiplicities[loop_id]
            push!(loop_data, (loop.weight, mult))
        end
        sort!(loop_data)
        
        # Check total weight matches cluster weight
        total_weight = sum(weight * mult for (weight, mult) in loop_data)
        if total_weight != cluster.weight
            push!(structure_violations, (site, cluster, "Weight mismatch: calculated=$total_weight, stored=$(cluster.weight)"))
            continue
        end
        
        # Check specific constraints based on weight
        if cluster.weight == 4
            # Weight-4 clusters should be single loops with weight 4, multiplicity 1
            if length(loop_data) != 1 || loop_data[1] != (4, 1)
                push!(structure_violations, (site, cluster, "Invalid weight-4 structure: $loop_data"))
            end
        elseif cluster.weight == 6
            # Weight-6 clusters should be single loops with weight 6, multiplicity 1
            if length(loop_data) != 1 || loop_data[1] != (6, 1)
                push!(structure_violations, (site, cluster, "Invalid weight-6 structure: $loop_data"))
            end
        end
    end
    
    if isempty(structure_violations)
        println("✅ All cluster structures valid: PASS")
        structure_ok = true
    else
        println("❌ Found $(length(structure_violations)) structure violations:")
        for (site, cluster, error_msg) in structure_violations
            println("  Site $site: $error_msg")
        end
        structure_ok = false
    end
    
    # Summary
    println("\n📊 MULTIPLICITY TEST SUMMARY:")
    println("="^30)
    println("✅ Multiplicity ≤ 2: $(multiplicity_ok ? "PASS" : "FAIL")")
    println("✅ Cluster structures valid: $(structure_ok ? "PASS" : "FAIL")")
    
    overall_pass = multiplicity_ok && structure_ok
    println("🏁 Overall result: $(overall_pass ? "✅ PASS" : "❌ FAIL")")
    
    return overall_pass
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_multiplicity_constraints()
end