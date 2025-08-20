#!/usr/bin/env julia

"""
debug_canonical_loops.jl

Debug why canonical loop representation isn't working for deduplication.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Serialization

function load_latest_11x11_file()
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    matching_files = filter(f -> contains(f, "L11") && contains(f, "w4") && contains(f, "periodic") && endswith(f, ".jld2"), files)
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end

function debug_canonical_representations()
    println("🔍 Debugging Canonical Loop Representations")
    println("="^50)
    
    data = load_latest_11x11_file()
    
    # Get first few weight-4 clusters from site 5
    site5_clusters = []
    for cluster in data.clusters_by_site[5]
        if cluster.weight == 4
            push!(site5_clusters, cluster)
        end
    end
    
    println("📊 Found $(length(site5_clusters)) weight-4 clusters on site 5")
    
    # Look at the loops in detail
    println("\n🔍 Loop details:")
    for (i, cluster) in enumerate(site5_clusters[1:min(4, end)])
        loop_id = cluster.loop_ids[1]
        loop = data.all_loops[loop_id]
        canonical = canonical_loop_representation(loop)
        
        println("  Cluster $i:")
        println("    Loop ID: $loop_id")
        println("    Vertices: $(loop.vertices)")
        println("    Edges: $(loop.edges)")
        println("    Canonical: $canonical")
        println()
    end
    
    # Check if these are actually the same loop translated
    if length(site5_clusters) >= 2
        loop1 = data.all_loops[site5_clusters[1].loop_ids[1]]
        loop2 = data.all_loops[site5_clusters[2].loop_ids[1]]
        
        println("🔍 Comparing first two loops:")
        println("  Loop 1 vertices: $(loop1.vertices)")
        println("  Loop 2 vertices: $(loop2.vertices)")
        
        # Check if they're translates of each other
        v1_sorted = sort(loop1.vertices)
        v2_sorted = sort(loop2.vertices)
        
        # Convert to coordinates to see the pattern
        L = 11
        function site_to_coords(site)
            i = div(site - 1, L) + 1
            j = mod(site - 1, L) + 1
            return (i, j)
        end
        
        coords1 = [site_to_coords(v) for v in v1_sorted]
        coords2 = [site_to_coords(v) for v in v2_sorted]
        
        println("  Loop 1 coords: $coords1")
        println("  Loop 2 coords: $coords2")
        
        # Check coordinate differences
        if length(coords1) == length(coords2)
            diffs = []
            for k in 1:length(coords1)
                di = coords2[k][1] - coords1[k][1]
                dj = coords2[k][2] - coords1[k][2]
                push!(diffs, (di, dj))
            end
            println("  Coordinate differences: $diffs")
            
            # Check if all differences are the same (translation)
            if length(Set(diffs)) == 1
                println("  ✅ These are translates of each other!")
                diff = diffs[1]
                println("  Translation: $diff")
            else
                println("  ❌ These are NOT simple translates")
            end
        end
    end
    
    # Test the canonical_loop_representation function directly
    println("\n🧪 Testing canonical_loop_representation:")
    
    # Create a simple test case
    if length(site5_clusters) >= 2
        loop1 = data.all_loops[site5_clusters[1].loop_ids[1]]
        loop2 = data.all_loops[site5_clusters[2].loop_ids[1]]
        
        canon1 = canonical_loop_representation(loop1)
        canon2 = canonical_loop_representation(loop2)
        
        println("  Loop 1 canonical: $canon1")
        println("  Loop 2 canonical: $canon2")
        println("  Are they equal? $(canon1 == canon2)")
        
        if canon1 != canon2
            println("  ❌ Canonical representations differ - this is the bug!")
            
            # Let's understand WHY they differ
            vertices1, edges1, weight1 = canon1
            vertices2, edges2, weight2 = canon2
            
            println("    Canon 1: vertices=$vertices1, edges=$edges1")
            println("    Canon 2: vertices=$vertices2, edges=$edges2")
        end
    end
    
    # Count how many unique canonical representations we have
    println("\n📊 Counting unique canonical representations:")
    canonical_set = Set()
    
    # Get all weight-4 loops
    weight4_loops = []
    for loop in data.all_loops
        if loop.weight == 4
            push!(weight4_loops, loop)
        end
    end
    
    println("  Total weight-4 loops: $(length(weight4_loops))")
    
    for loop in weight4_loops
        canonical = canonical_loop_representation(loop)
        push!(canonical_set, canonical)
    end
    
    println("  Unique canonical forms: $(length(canonical_set))")
    
    if length(canonical_set) == 1
        println("  ✅ All weight-4 loops have the same canonical form!")
    elseif length(canonical_set) == 4
        println("  🤔 4 unique forms - might be different loop types (2x2 squares in different orientations?)")
    else
        println("  ❌ Unexpected number of canonical forms")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    debug_canonical_representations()
end