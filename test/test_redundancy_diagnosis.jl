#!/usr/bin/env julia

"""
test_redundancy_diagnosis.jl

Diagnose why there are 484 weight-4 multiplicity-2 clusters instead of 121.
"""

# Include necessary modules
include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Serialization

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

function diagnose_weight4_redundancy()
    """Diagnose the weight-4 multiplicity-2 cluster redundancy issue."""
    println("🔍 Diagnosing Weight-4 Cluster Redundancy Issue")
    println("="^50)
    
    # Load data
    data = load_latest_cluster_file()
    L = data.lattice_size
    
    println("📊 Lattice: $(L)×$(L), Expected unique clusters: $(L*L) = 121")
    
    # Find all multiplicity-2 clusters (single loop, mult=2)
    mult2_clusters = []
    
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            # Check if this has any loop with multiplicity 2
            for (loop_id, mult) in cluster.multiplicities
                if mult == 2
                    loop = data.all_loops[loop_id]
                    push!(mult2_clusters, (site, cluster, loop_id, loop))
                end
            end
        end
    end
    
    println("📊 Found $(length(mult2_clusters)) multiplicity-2 clusters")
    
    if length(mult2_clusters) != 484
        println("⚠️  Expected 484 based on test result, got $(length(mult2_clusters))")
    end
    
    # Analyze the underlying loops
    unique_loops = Set{Int}()
    loop_usage_count = Dict{Int, Int}()
    loop_weights = Dict{Int, Int}()
    
    for (site, cluster, loop_id, loop) in mult2_clusters
        push!(unique_loops, loop_id)
        loop_usage_count[loop_id] = get(loop_usage_count, loop_id, 0) + 1
        loop_weights[loop_id] = loop.weight
    end
    
    println("📊 These clusters use $(length(unique_loops)) unique loops")
    println("📊 Expected: 4 unique loops (since 484 = 4 × 121)")
    
    # Show the usage pattern
    println("\n🔍 Loop usage pattern:")
    for loop_id in sort(collect(unique_loops))[1:min(10, length(unique_loops))]  # Show first 10
        count = loop_usage_count[loop_id]
        weight = loop_weights[loop_id]
        loop = data.all_loops[loop_id]
        println("  Loop $loop_id (weight $weight): used in $count clusters")
        println("    Vertices: $(loop.vertices)")
        println("    Edges: $(loop.edges)")
    end
    
    if length(unique_loops) > 10
        println("  ... (showing first 10 of $(length(unique_loops)) loops)")
    end
    
    # Check if there should be only 1 unique loop
    if length(unique_loops) > 1
        println("\n🔍 Analyzing if these are actually the same loop under translation...")
        
        # Take the first loop as reference
        first_loop_id = first(unique_loops)
        reference_loop = data.all_loops[first_loop_id]
        
        println("  Reference loop (ID $first_loop_id): $(reference_loop.vertices)")
        
        # Check canonical representations
        println("\n🔍 Canonical representations:")
        for loop_id in sort(collect(unique_loops))
            loop = data.all_loops[loop_id]
            canonical = canonical_loop_representation(loop)
            println("  Loop $loop_id: $canonical")
        end
        
        # Check if they have the same canonical representation
        canonical_forms = Set{Tuple}()
        for loop_id in unique_loops
            loop = data.all_loops[loop_id]
            canonical = canonical_loop_representation(loop)
            push!(canonical_forms, canonical)
        end
        
        println("📊 Unique canonical forms: $(length(canonical_forms))")
        if length(canonical_forms) == 1
            println("✅ All loops have the same canonical form - redundancy issue confirmed!")
        else
            println("❌ Different canonical forms found - these are genuinely different loops")
        end
    end
    
    # Analyze cluster distribution by site
    println("\n🔍 Analyzing cluster distribution by site...")
    
    site_counts = Dict{Int, Int}()
    for (site, cluster, loop_id, loop) in mult2_clusters
        site_counts[site] = get(site_counts, site, 0) + 1
    end
    
    counts_histogram = Dict{Int, Int}()
    for count in values(site_counts)
        counts_histogram[count] = get(counts_histogram, count, 0) + 1
    end
    
    println("📊 Clusters per site histogram:")
    for (count, num_sites) in sort(collect(counts_histogram))
        println("  $count clusters: $num_sites sites")
    end
    
    expected_count_per_site = 484 ÷ 121  # Should be 4 if evenly distributed
    println("📊 Expected clusters per site: $expected_count_per_site")
    
    # Sample detailed analysis
    println("\n🔍 Sample detailed analysis (first 5 clusters):")
    for i in 1:min(5, length(mult2_clusters))
        site, cluster, loop_id, loop = mult2_clusters[i]
        println("  $i. Site $site:")
        println("     Loop ID: $loop_id")
        println("     Loop weight: $(loop.weight)")
        println("     Loop vertices: $(loop.vertices)")
        println("     Cluster signature: $(canonical_cluster_signature(cluster))")
    end
    
    return length(mult2_clusters), length(unique_loops)
end

# Run the diagnosis
if abspath(PROGRAM_FILE) == @__FILE__
    diagnose_weight4_redundancy()
end