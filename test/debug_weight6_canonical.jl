#!/usr/bin/env julia

"""
debug_weight6_canonical.jl

Debug the canonical form generation for weight-6 loops to understand why deduplication 
is not working properly (442 clusters instead of expected 242).
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Serialization

function load_latest_w6_file()
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    matching_files = filter(f -> contains(f, "L11") && contains(f, "w6") && contains(f, "test_w6") && endswith(f, ".jld2"), files)
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end

function debug_weight6_canonical_forms()
    println("🔍 Debugging Weight-6 Canonical Forms")
    println("="^50)
    
    data = load_latest_w6_file()
    L = data.lattice_size
    
    # Get all weight-6 loops
    weight6_loops = [loop for loop in data.all_loops if loop.weight == 6]
    println("📊 Total weight-6 loops: $(length(weight6_loops))")
    
    # Test current canonical form function
    println("\n🧪 Testing current canonical form generation...")
    
    # Include the coordinate-based canonical form function from generate_ising_clusters.jl
    function create_coordinate_canonical_form(loop, L::Int = 11)
        # Convert vertices to coordinates
        function site_to_coords(site::Int)
            i = div(site - 1, L) + 1
            j = mod(site - 1, L) + 1
            return (i, j)
        end
        
        coords = [site_to_coords(v) for v in sort(loop.vertices)]
        
        # Normalize coordinates to start from (0, 0)
        min_i = minimum(coord[1] for coord in coords)
        min_j = minimum(coord[2] for coord in coords)
        
        normalized_coords = [(coord[1] - min_i, coord[2] - min_j) for coord in coords]
        
        # Handle periodic wrapping for 2x2 squares (all weight-4 loops)
        if length(normalized_coords) == 4 && loop.weight == 4
            i_values = [coord[1] for coord in coords]
            j_values = [coord[2] for coord in coords]
            
            i_span = maximum(i_values) - minimum(i_values)
            j_span = maximum(j_values) - minimum(j_values)
            
            # If span > L/2, it's likely wrapped around due to periodic boundaries
            if i_span > L/2 || j_span > L/2
                # All wrapped 2x2 squares have the same canonical form
                normalized_coords = [(0, 0), (0, 1), (1, 0), (1, 1)]
            end
        end
        
        # Sort for canonical ordering
        normalized_coords = sort(normalized_coords)
        
        return (normalized_coords, loop.weight)
    end
    
    # Test canonical forms on weight-6 loops
    canonical_forms = Dict()
    
    for (i, loop) in enumerate(weight6_loops[1:min(10, end)])  # Test first 10
        canonical = create_coordinate_canonical_form(loop, L)
        
        if !haskey(canonical_forms, canonical)
            canonical_forms[canonical] = []
        end
        push!(canonical_forms[canonical], (i, loop))
        
        println("  Loop $i: vertices=$(loop.vertices)")
        
        # Convert to coordinates
        function site_to_coords(site::Int)
            i = div(site - 1, L) + 1
            j = mod(site - 1, L) + 1
            return (i, j)
        end
        
        coords = [site_to_coords(v) for v in sort(loop.vertices)]
        println("    Coordinates: $coords")
        println("    Canonical form: $canonical")
        println()
    end
    
    println("📊 Unique canonical forms in sample: $(length(canonical_forms))")
    
    # Test all weight-6 loops
    println("\n🔄 Testing all $(length(weight6_loops)) weight-6 loops...")
    all_canonical_forms = Set()
    
    for loop in weight6_loops
        canonical = create_coordinate_canonical_form(loop, L)
        push!(all_canonical_forms, canonical)
    end
    
    println("📊 Total unique canonical forms: $(length(all_canonical_forms))")
    
    if length(all_canonical_forms) == 2
        println("  ✅ Perfect! This matches the expected 2 unique weight-6 loop types")
        println("  Expected clusters: 2 types × 121 sites = 242 clusters")
    elseif length(all_canonical_forms) > 2
        println("  🤔 More than 2 unique forms. Let's investigate the different types:")
        
        # Group loops by canonical form
        loops_by_canonical = Dict()
        for loop in weight6_loops
            canonical = create_coordinate_canonical_form(loop, L)
            if !haskey(loops_by_canonical, canonical)
                loops_by_canonical[canonical] = []
            end
            push!(loops_by_canonical[canonical], loop)
        end
        
        for (canonical, loops) in loops_by_canonical
            println("    Form $canonical: $(length(loops)) loops")
            if length(loops) > 0
                example_loop = loops[1]
                coords = [site_to_coords(v) for v in sort(example_loop.vertices)]
                println("      Example coordinates: $coords")
            end
        end
        
        println("  🔍 This explains why we get more clusters than expected!")
    end
    
    # Let's also check what types of geometric patterns exist for weight-6
    println("\n🔍 Analyzing weight-6 loop geometric patterns:")
    
    function analyze_loop_pattern(loop, L)
        function site_to_coords(site::Int)
            i = div(site - 1, L) + 1
            j = mod(site - 1, L) + 1
            return (i, j)
        end
        
        coords = [site_to_coords(v) for v in sort(loop.vertices)]
        
        # Calculate spans
        i_values = [c[1] for c in coords]
        j_values = [c[2] for c in coords]
        
        i_span = maximum(i_values) - minimum(i_values)
        j_span = maximum(j_values) - minimum(j_values)
        
        # Check if it's a rectangular pattern
        is_rect = length(Set(i_values)) <= 3 && length(Set(j_values)) <= 3
        
        return (length(coords), i_span, j_span, is_rect)
    end
    
    pattern_types = Dict()
    
    for loop in weight6_loops[1:min(50, end)]  # Analyze first 50
        pattern = analyze_loop_pattern(loop, L)
        pattern_types[pattern] = get(pattern_types, pattern, 0) + 1
    end
    
    println("  Pattern analysis (first 50 loops):")
    for (pattern, count) in sort(collect(pattern_types), by=x->x[2], rev=true)
        n_vertices, i_span, j_span, is_rect = pattern
        println("    $count loops: $n_vertices vertices, spans=($i_span,$j_span), rectangular=$is_rect")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    debug_weight6_canonical_forms()
end