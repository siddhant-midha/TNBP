#!/usr/bin/env julia

"""
check_vertex_uniqueness.jl

Check if each vertex has a unique index and understand the loop structure better.
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

function check_vertex_uniqueness()
    println("🔍 Checking Vertex Uniqueness and Loop Structure")
    println("="^60)
    
    data = load_latest_11x11_file()
    L = data.lattice_size
    println("📊 Lattice size: $L×$L = $(L*L) vertices")
    
    # Get weight-4 clusters from site 5 to examine
    site5_clusters = []
    for cluster in data.clusters_by_site[5]
        if cluster.weight == 4
            push!(site5_clusters, cluster)
        end
    end
    
    println("📊 Found $(length(site5_clusters)) weight-4 clusters on site 5")
    
    # Convert site number to coordinates
    function site_to_coords(site::Int)
        i = div(site - 1, L) + 1
        j = mod(site - 1, L) + 1
        return (i, j)
    end
    
    function coords_to_site(i::Int, j::Int)
        # Handle periodic wrapping
        i_wrap = mod(i - 1, L) + 1
        j_wrap = mod(j - 1, L) + 1
        return (i_wrap - 1) * L + j_wrap
    end
    
    println("\n🔍 Detailed analysis of loops:")
    
    for (idx, cluster) in enumerate(site5_clusters)
        loop_id = cluster.loop_ids[1]
        loop = data.all_loops[loop_id]
        
        println("  Loop $idx (ID: $loop_id):")
        println("    Vertices: $(loop.vertices)")
        
        # Convert to coordinates
        coords = [site_to_coords(v) for v in loop.vertices]
        println("    Coordinates: $coords")
        
        # Check if vertices are unique
        unique_vertices = Set(loop.vertices)
        println("    Unique vertices: $(length(unique_vertices) == length(loop.vertices) ? "✅ YES" : "❌ NO")")
        
        # Check if this forms a connected 2×2 square
        i_coords = [c[1] for c in coords]
        j_coords = [c[2] for c in coords]
        
        println("    i-coordinates: $i_coords (range: $(minimum(i_coords))-$(maximum(i_coords)))")
        println("    j-coordinates: $j_coords (range: $(minimum(j_coords))-$(maximum(j_coords)))")
        
        # Check edges
        println("    Edges: $(loop.edges)")
        
        # Verify edges connect the vertices properly
        edge_vertices = Set()
        for (u, v) in loop.edges
            push!(edge_vertices, u)
            push!(edge_vertices, v)
        end
        vertices_from_edges = Set(loop.vertices)
        
        if edge_vertices == vertices_from_edges
            println("    Edge consistency: ✅ All vertices appear in edges")
        else
            println("    Edge consistency: ❌ Mismatch between vertices and edges")
            println("      Vertices: $vertices_from_edges")
            println("      From edges: $edge_vertices")
        end
        
        println()
    end
    
    # Now let's check: are these actually different loops or the same loop?
    println("🤔 Analysis: Are these different loops?")
    
    if length(site5_clusters) >= 2
        loop1 = data.all_loops[site5_clusters[1].loop_ids[1]]
        loop3 = data.all_loops[site5_clusters[3].loop_ids[1]]
        
        coords1 = [site_to_coords(v) for v in sort(loop1.vertices)]
        coords3 = [site_to_coords(v) for v in sort(loop3.vertices)]
        
        println("  Loop 1 coordinates (sorted): $coords1")
        println("  Loop 3 coordinates (sorted): $coords3")
        
        # Check if they form the same geometric pattern
        # For 2x2 squares, normalize to relative coordinates
        function normalize_coords(coords)
            min_i = minimum(c[1] for c in coords)
            min_j = minimum(c[2] for c in coords) 
            return [(c[1] - min_i, c[2] - min_j) for c in coords]
        end
        
        norm1 = sort(normalize_coords(coords1))
        norm3 = sort(normalize_coords(coords3))
        
        println("  Loop 1 normalized: $norm1")
        println("  Loop 3 normalized: $norm3")
        
        if norm1 == norm3
            println("  ✅ Same geometric pattern - these should be equivalent!")
        else
            println("  ❌ Different geometric patterns - these are truly different loops")
        end
    end
    
    # Check all weight-4 loops to see how many unique geometric patterns exist
    println("\n📊 Checking all weight-4 loops for unique patterns:")
    
    weight4_loops = [loop for loop in data.all_loops if loop.weight == 4]
    println("  Total weight-4 loops: $(length(weight4_loops))")
    
    # Group by geometric pattern
    pattern_groups = Dict()
    
    for loop in weight4_loops
        coords = [site_to_coords(v) for v in sort(loop.vertices)]
        
        # Normalize coordinates
        min_i = minimum(c[1] for c in coords)
        min_j = minimum(c[2] for c in coords)
        normalized = sort([(c[1] - min_i, c[2] - min_j) for c in coords])
        
        # Handle periodic wrapping for 2x2 squares
        if length(normalized) == 4
            i_span = maximum(c[1] for c in coords) - minimum(c[1] for c in coords)
            j_span = maximum(c[2] for c in coords) - minimum(c[2] for c in coords)
            
            # If span > L/2, it's likely wrapped around
            if i_span > L/2 || j_span > L/2
                # Normalize wrapped squares to standard 2x2
                normalized = [(0, 0), (0, 1), (1, 0), (1, 1)]
            end
        end
        
        pattern_key = tuple(normalized...)
        if !haskey(pattern_groups, pattern_key)
            pattern_groups[pattern_key] = []
        end
        push!(pattern_groups[pattern_key], loop)
    end
    
    println("  Unique geometric patterns: $(length(pattern_groups))")
    
    for (pattern, loops) in pattern_groups
        println("    Pattern $pattern: $(length(loops)) loops")
    end
    
    if length(pattern_groups) == 1
        println("  🎯 All weight-4 loops have the same geometric pattern!")
        println("  This confirms they should be deduplicated to 121 clusters")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    check_vertex_uniqueness()
end