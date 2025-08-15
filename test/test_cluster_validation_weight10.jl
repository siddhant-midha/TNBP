#!/usr/bin/env julia

"""
test_cluster_validation_weight10.jl

Comprehensive validation test for cluster enumeration results.
Tests weight 10, size 11, PBC cluster enumeration for:
1. Multiplicity constraints (≤ 2)
2. Loop weight constraints for multiplicity=2 clusters ((4,4) or (4,6))
3. Vertex sharing constraints (1, 2, or 4 vertices shared)
4. Translation symmetry (each equivalence class has exactly 11² = 121 clusters)
"""

# Include necessary modules
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("../functions/loop_enumeration_square.jl")

using Test
using ProgressMeter
using Serialization

# Translation utilities (from generate_ising_clusters_one_site_old.jl)
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

function count_shared_vertices(loop1::Loop, loop2::Loop)
    """Count the number of shared vertices between two loops."""
    vertices1 = Set(loop1.vertices)
    vertices2 = Set(loop2.vertices)
    return length(intersect(vertices1, vertices2))
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

function validate_multiplicity_constraints(data::ClusterEnumerationData)
    """Test 1: Check that all clusters have multiplicity ≤ 2."""
    println("\n🔍 Test 1: Validating multiplicity constraints...")
    
    violation_count = 0
    total_clusters = 0
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            total_clusters += 1
            max_multiplicity = maximum(values(cluster.multiplicities))
            
            if max_multiplicity > 2
                violation_count += 1
                println("  ❌ VIOLATION: Site $site has cluster with multiplicity $max_multiplicity")
                println("     Cluster: $(cluster.loop_ids) with multiplicities $(cluster.multiplicities)")
                
                if violation_count >= 5  # Limit output
                    println("     ... (limiting output to first 5 violations)")
                    break
                end
            end
        end
        
        if violation_count >= 5
            break
        end
    end
    
    println("  📊 Checked $total_clusters clusters")
    if violation_count == 0
        println("  ✅ All clusters have multiplicity ≤ 2")
    else
        println("  ❌ Found $violation_count violations")
    end
    
    return violation_count == 0
end

function validate_loop_weight_constraints(data::ClusterEnumerationData)
    """Test 2: For multiplicity=2 clusters, check loop weights are (4,4) or (4,6)."""
    println("\n🔍 Test 2: Validating loop weight constraints for multiplicity=2 clusters...")
    
    violation_count = 0
    mult2_count = 0
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            # Check if this cluster has any loop with multiplicity 2
            has_mult2 = any(mult == 2 for mult in values(cluster.multiplicities))
            
            if has_mult2
                mult2_count += 1
                
                # Get the loops with multiplicity 2
                mult2_loops = [loop_id for (loop_id, mult) in cluster.multiplicities if mult == 2]
                
                for loop_id in mult2_loops
                    loop_weight = data.all_loops[loop_id].weight
                    
                    # Check if this is a valid two-loop cluster
                    if length(cluster.loop_ids) == 1
                        # Single loop with multiplicity 2 - check weight is 4
                        if loop_weight != 4
                            violation_count += 1
                            println("  ❌ VIOLATION: Single loop (mult=2) has weight $loop_weight, expected 4")
                            
                            if violation_count >= 5
                                break
                            end
                        end
                    elseif length(cluster.loop_ids) == 2
                        # Two different loops - check weights are (4,4) or (4,6)
                        weights = [data.all_loops[id].weight for id in cluster.loop_ids]
                        sort!(weights)
                        
                        if !(weights == [4, 4] || weights == [4, 6])
                            violation_count += 1
                            println("  ❌ VIOLATION: Two-loop cluster has weights $weights, expected [4,4] or [4,6]")
                            
                            if violation_count >= 5
                                break
                            end
                        end
                    end
                end
                
                if violation_count >= 5
                    break
                end
            end
        end
        
        if violation_count >= 5
            break
        end
    end
    
    println("  📊 Checked $mult2_count clusters with multiplicity=2")
    if violation_count == 0
        println("  ✅ All multiplicity=2 clusters have valid loop weights")
    else
        println("  ❌ Found $violation_count violations")
    end
    
    return violation_count == 0
end

function validate_vertex_sharing_constraints(data::ClusterEnumerationData)
    """Test 3: For multiplicity=2 clusters with two different loops, check they share 1, 2, or 4 vertices."""
    println("\n🔍 Test 3: Validating vertex sharing constraints...")
    
    violation_count = 0
    checked_count = 0
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            # Only check clusters with two different loops
            if length(cluster.loop_ids) == 2
                checked_count += 1
                
                loop1 = data.all_loops[cluster.loop_ids[1]]
                loop2 = data.all_loops[cluster.loop_ids[2]]
                
                shared_vertices = count_shared_vertices(loop1, loop2)
                
                if !(shared_vertices in [1, 2, 4])
                    violation_count += 1
                    println("  ❌ VIOLATION: Two loops share $shared_vertices vertices, expected 1, 2, or 4")
                    println("     Loop 1: $(loop1.vertices)")
                    println("     Loop 2: $(loop2.vertices)")
                    
                    if violation_count >= 5
                        break
                    end
                end
            end
        end
        
        if violation_count >= 5
            break
        end
    end
    
    println("  📊 Checked $checked_count two-loop clusters")
    if violation_count == 0
        println("  ✅ All two-loop clusters have valid vertex sharing (1, 2, or 4 vertices)")
    else
        println("  ❌ Found $violation_count violations")
    end
    
    return violation_count == 0
end

function validate_translation_symmetry(data::ClusterEnumerationData)
    """Test 4: Check that each translation equivalence class has exactly L² = 121 clusters."""
    println("\n🔍 Test 4: Validating translation symmetry...")
    
    L = data.lattice_size
    expected_class_size = L * L  # 11² = 121
    
    # Collect all clusters from all sites
    all_clusters = Cluster[]
    for clusters in values(data.clusters_by_site)
        append!(all_clusters, clusters)
    end
    
    println("  📊 Processing $(length(all_clusters)) total clusters...")
    
    # Group clusters by their canonical representative under translation
    equivalence_classes = Dict{Tuple, Vector{Cluster}}()
    seen_signatures = Set{Tuple}()
    
    progress = Progress(length(all_clusters), dt=1.0, desc="Grouping by translation: ", color=:blue)
    
    for cluster in all_clusters
        # Find canonical representative
        canonical_cluster = canonical_cluster_under_translation(cluster, L, data.all_loops)
        canonical_signature = canonical_cluster_signature(canonical_cluster)
        
        # Add to equivalence class
        if !haskey(equivalence_classes, canonical_signature)
            equivalence_classes[canonical_signature] = Cluster[]
        end
        push!(equivalence_classes[canonical_signature], cluster)
        
        next!(progress)
    end
    
    finish!(progress)
    
    # Check each equivalence class size
    violation_count = 0
    class_sizes = Int[]
    
    println("  📊 Found $(length(equivalence_classes)) equivalence classes")
    
    for (signature, class_clusters) in equivalence_classes
        class_size = length(class_clusters)
        push!(class_sizes, class_size)
        
        if class_size != expected_class_size
            violation_count += 1
            
            if violation_count <= 5  # Limit output
                println("  ❌ VIOLATION: Equivalence class has $class_size clusters, expected $expected_class_size")
                println("     Canonical signature: $signature")
            end
        end
    end
    
    # Statistics
    min_size = minimum(class_sizes)
    max_size = maximum(class_sizes)
    unique_sizes = sort(unique(class_sizes))
    
    println("  📊 Class size statistics:")
    println("    Min: $min_size, Max: $max_size")
    println("    Unique sizes: $unique_sizes")
    
    if violation_count == 0
        println("  ✅ All equivalence classes have exactly $expected_class_size clusters")
    else
        println("  ❌ Found $violation_count classes with incorrect size")
    end
    
    return violation_count == 0
end

function run_comprehensive_validation()
    """Run all validation tests."""
    println("🧪 Comprehensive Cluster Validation Test")
    println("="^60)
    println("Testing weight 10, size 11, PBC cluster enumeration results")
    
    # Load data
    local data
    try
        data = load_latest_cluster_file()
        println("✅ Successfully loaded data")
        
        # Print summary
        total_clusters = sum(length(clusters) for clusters in values(data.clusters_by_site))
        println("📊 Summary: $(data.lattice_size)×$(data.lattice_size) lattice, max weight $(data.max_weight)")
        println("📊 Total clusters: $total_clusters, Total loops: $(length(data.all_loops))")
        
    catch e
        println("❌ Failed to load cluster data: $e")
        return false
    end
    
    # Run all tests
    results = Bool[]
    
    push!(results, validate_multiplicity_constraints(data))
    push!(results, validate_loop_weight_constraints(data))
    push!(results, validate_vertex_sharing_constraints(data))
    push!(results, validate_translation_symmetry(data))
    
    # Summary
    println("\n" * "="^60)
    println("🏁 VALIDATION SUMMARY")
    println("="^60)
    
    test_names = [
        "Multiplicity constraints (≤ 2)",
        "Loop weight constraints ((4,4) or (4,6))", 
        "Vertex sharing constraints (1, 2, or 4 vertices)",
        "Translation symmetry (121 clusters per class)"
    ]
    
    all_passed = true
    for (i, (name, passed)) in enumerate(zip(test_names, results))
        status = passed ? "✅ PASS" : "❌ FAIL"
        println("Test $i: $status - $name")
        all_passed = all_passed && passed
    end
    
    println("\n" * "="^60)
    if all_passed
        println("🎉 ALL TESTS PASSED! Cluster enumeration is valid.")
    else
        println("⚠️  SOME TESTS FAILED! Review the violations above.")
    end
    
    return all_passed
end

# Run the test if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    @testset "Cluster Validation Weight 10" begin
        @test run_comprehensive_validation()
    end
end