"""
LoopEnumeration.jl

Julia implementation for enumerating loops (connected subgraphs with degree ≥ 2)
on periodic square lattices.
"""

struct Loop
    vertices::Vector{Int}
    edges::Vector{Tuple{Int,Int}}
    weight::Int
end

# Canonical representation for loop deduplication
function canonical_loop_representation(loop::Loop)
    vertices = sort(loop.vertices)  # Keep ALL vertices
    edges = sort([(min(u, v), max(u, v)) for (u, v) in loop.edges])  # Convert to tuples properly
    return (vertices, edges, loop.weight)
end

struct LoopEnumerator
    adj_matrix::Matrix{Int}
    n_vertices::Int
    adj_list::Dict{Int, Vector{Int}}
    L::Int  # Lattice size for locality optimization
    coords::Dict{Int, Tuple{Int, Int}}  # site -> (i, j) coordinates
    
    function LoopEnumerator(adj_matrix::Matrix{Int})
        n = size(adj_matrix, 1)
        L = Int(sqrt(n))  # Assume square lattice
        adj_list = build_adjacency_list(adj_matrix, n)
        coords = build_coordinate_map(L)
        new(adj_matrix, n, adj_list, L, coords)
    end
end

function build_adjacency_list(adj_matrix::Matrix{Int}, n::Int)
    adj_list = Dict{Int, Vector{Int}}()
    for i in 1:n
        adj_list[i] = Int[]
        for j in 1:n
            if adj_matrix[i, j] != 0
                push!(adj_list[i], j)
            end
        end
    end
    return adj_list
end

function build_coordinate_map(L::Int)
    """Build mapping from site index to (i,j) coordinates."""
    coords = Dict{Int, Tuple{Int, Int}}()
    for i in 1:L, j in 1:L
        site = (i-1) * L + j
        coords[site] = (i, j)
    end
    return coords
end

function lattice_distance(coord1::Tuple{Int,Int}, coord2::Tuple{Int,Int}, L::Int)
    """Compute shortest distance on periodic square lattice."""
    i1, j1 = coord1
    i2, j2 = coord2
    
    di = min(abs(i2 - i1), L - abs(i2 - i1))
    dj = min(abs(j2 - j1), L - abs(j2 - j1))
    
    return di + dj  # Manhattan distance with periodic boundaries
end

function find_loops_supported_on_vertex(enumerator::LoopEnumerator, vertex::Int, max_weight::Int)
    """Find all loops supported on a vertex using DFS."""
    found_loops = Loop[]
    
    dfs_loop_enumeration!(enumerator, vertex, max_weight, found_loops)
    return found_loops
end

function find_loops_supported_on_vertex_optimized(enumerator::LoopEnumerator, vertex::Int, max_weight::Int, 
                                                  global_seen::Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}})
    """Find loops supported on vertex with minimal but effective optimization."""
    found_loops = Loop[]
    
    # Use the lightweight optimized version
    dfs_loop_enumeration_optimized!(enumerator, vertex, max_weight, found_loops, global_seen)
    return found_loops
end

function find_loops_supported_on_vertex_fast(enumerator::LoopEnumerator, vertex::Int, max_weight::Int)
    """Ultra-minimal optimization - just canonical filtering during enumeration."""
    found_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    # Get all edges adjacent to vertex
    center_edges = Tuple{Int,Int}[]
    for neighbor in enumerator.adj_list[vertex]
        edge = (min(vertex, neighbor), max(vertex, neighbor))
        push!(center_edges, edge)
    end
    
    function dfs_build_edge_set(current_edges::Vector{Tuple{Int,Int}}, edge_index::Int)
        # Check if current edge set forms a valid loop
        if length(current_edges) >= 2
            vertices_set = Set{Int}()
            for (u, v) in current_edges
                push!(vertices_set, u)
                push!(vertices_set, v)
            end
            
            if (vertex in vertices_set && is_valid_loop(vertices_set, current_edges))
                # ONLY optimization: canonical filtering
                if vertex == minimum(vertices_set)
                    loop = Loop(sort(collect(vertices_set)), sort(current_edges), length(current_edges))
                    canonical = canonical_loop_representation(loop)
                    
                    if !(canonical in seen_loops)
                        push!(seen_loops, canonical)
                        push!(found_loops, loop)
                    end
                end
            end
        end
        
        # Simple early stopping
        if length(current_edges) >= max_weight
            return
        end
        
        # Continue search
        current_vertices = Set{Int}()
        for (u, v) in current_edges
            push!(current_vertices, u)
            push!(current_vertices, v)
        end
        
        if isempty(current_vertices)
            for edge in center_edges
                dfs_build_edge_set([edge], 1)
            end
        else
            candidate_edges = Set{Tuple{Int,Int}}()
            for v in current_vertices
                for neighbor in enumerator.adj_list[v]
                    edge = (min(v, neighbor), max(v, neighbor))
                    if !(edge in current_edges)
                        push!(candidate_edges, edge)
                    end
                end
            end
            
            for edge in candidate_edges
                new_edges = copy(current_edges)
                push!(new_edges, edge)
                dfs_build_edge_set(new_edges, edge_index + 1)
            end
        end
    end
    
    dfs_build_edge_set(Tuple{Int,Int}[], 0)
    return found_loops
end

function dfs_loop_enumeration!(enumerator::LoopEnumerator, center_vertex::Int, 
                              max_weight::Int, found_loops::Vector{Loop})
    """DFS-based loop enumeration: systematically explore all edge combinations."""
    
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    # Get all edges adjacent to center_vertex to start exploration
    center_edges = Tuple{Int,Int}[]
    for neighbor in enumerator.adj_list[center_vertex]
        edge = (min(center_vertex, neighbor), max(center_vertex, neighbor))
        push!(center_edges, edge)
    end
    
    # DFS to systematically build edge sets
    function dfs_build_edge_set(current_edges::Vector{Tuple{Int,Int}}, edge_index::Int)
        
        # Check if current edge set forms a valid loop
        if length(current_edges) >= 2
            # Extract vertices from edges
            vertices_set = Set{Int}()
            for (u, v) in current_edges
                push!(vertices_set, u)
                push!(vertices_set, v)
            end
            
            # Check if center_vertex is included and if it's a valid loop
            if (center_vertex in vertices_set && 
                is_valid_loop(vertices_set, current_edges))
                
                # Create loop and check for duplicates
                loop = Loop(sort(collect(vertices_set)), sort(current_edges), length(current_edges))
                canonical = canonical_loop_representation(loop)
                
                if !(canonical in seen_loops)
                    push!(seen_loops, canonical)
                    push!(found_loops, loop)
                end
            end
        end
        
        # Stop if we've reached max weight
        if length(current_edges) >= max_weight
            return
        end
        
        # Get all possible next edges from current vertices
        current_vertices = Set{Int}()
        for (u, v) in current_edges
            push!(current_vertices, u)
            push!(current_vertices, v)
        end
        
        if isempty(current_vertices)
            # Start with edges from center_vertex
            for edge in center_edges
                dfs_build_edge_set([edge], 1)
            end
        else
            # Add edges from current vertices
            candidate_edges = Set{Tuple{Int,Int}}()
            for vertex in current_vertices
                for neighbor in enumerator.adj_list[vertex]
                    edge = (min(vertex, neighbor), max(vertex, neighbor))
                    if !(edge in current_edges)  # Don't repeat edges
                        push!(candidate_edges, edge)
                    end
                end
            end
            
            # Try each candidate edge
            for edge in candidate_edges
                new_edges = copy(current_edges)
                push!(new_edges, edge)
                dfs_build_edge_set(new_edges, edge_index + 1)
            end
        end
    end
    
    # Start DFS with empty edge set
    dfs_build_edge_set(Tuple{Int,Int}[], 0)
end

function dfs_loop_enumeration_optimized!(enumerator::LoopEnumerator, center_vertex::Int, 
                                        max_weight::Int, found_loops::Vector{Loop},
                                        global_seen::Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}})
    """Lightweight optimized DFS - removes expensive operations causing slowdown."""
    
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    # Get all edges adjacent to center_vertex to start exploration
    center_edges = Tuple{Int,Int}[]
    for neighbor in enumerator.adj_list[center_vertex]
        edge = (min(center_vertex, neighbor), max(center_vertex, neighbor))
        push!(center_edges, edge)
    end
    
    # DFS to systematically build edge sets - similar to original but with key optimizations
    function dfs_build_edge_set(current_edges::Vector{Tuple{Int,Int}}, edge_index::Int)
        
        # Check if current edge set forms a valid loop
        if length(current_edges) >= 2
            # Extract vertices from edges
            vertices_set = Set{Int}()
            for (u, v) in current_edges
                push!(vertices_set, u)
                push!(vertices_set, v)
            end
            
            # Check if center_vertex is included and if it's a valid loop
            if (center_vertex in vertices_set && 
                is_valid_loop(vertices_set, current_edges))
                
                # OPTIMIZATION 1: Only proceed if center_vertex is minimum (reduces duplicates)
                if center_vertex == minimum(vertices_set)
                    # Create loop and check for duplicates
                    loop = Loop(sort(collect(vertices_set)), sort(current_edges), length(current_edges))
                    canonical = canonical_loop_representation(loop)
                    
                    # OPTIMIZATION 2: Skip expensive global_seen check for now - use local only
                    if !(canonical in seen_loops)
                        push!(seen_loops, canonical)
                        push!(found_loops, loop)
                    end
                end
            end
        end
        
        # OPTIMIZATION 3: Simple weight-based early stopping
        if length(current_edges) >= max_weight
            return
        end
        
        # Get all possible next edges from current vertices
        current_vertices = Set{Int}()
        for (u, v) in current_edges
            push!(current_vertices, u)
            push!(current_vertices, v)
        end
        
        if isempty(current_vertices)
            # Start with edges from center_vertex
            for edge in center_edges
                dfs_build_edge_set([edge], 1)
            end
        else
            # Add edges from current vertices - REMOVED expensive distance checks
            candidate_edges = Set{Tuple{Int,Int}}()
            for vertex in current_vertices
                for neighbor in enumerator.adj_list[vertex]
                    edge = (min(vertex, neighbor), max(vertex, neighbor))
                    if !(edge in current_edges)  # Don't repeat edges
                        push!(candidate_edges, edge)
                    end
                end
            end
            
            # Try each candidate edge
            for edge in candidate_edges
                new_edges = copy(current_edges)
                push!(new_edges, edge)
                dfs_build_edge_set(new_edges, edge_index + 1)
            end
        end
    end
    
    # Start DFS with empty edge set
    dfs_build_edge_set(Tuple{Int,Int}[], 0)
end

function could_be_valid_loop(vertices::Set{Int}, num_edges::Int)
    """Quick check if vertex count and edge count could form a valid loop."""
    n_vertices = length(vertices)
    
    # For a connected graph with all vertices having degree ≥ 2:
    # minimum edges needed is n_vertices (for a cycle)
    # maximum reasonable edges is roughly 2*n_vertices for dense loops
    return num_edges >= n_vertices && num_edges <= 3*n_vertices
end

function is_reasonable_extension(current_vertices::Set{Int}, new_vertex::Int, enumerator::LoopEnumerator)
    """Check if adding new_vertex creates a reasonable loop extension."""
    
    # For large loops, use spatial locality constraints
    if length(current_vertices) > 4
        # Get coordinates of current vertices and new vertex
        center_coords = enumerator.coords[minimum(current_vertices)]
        new_coords = enumerator.coords[new_vertex]
        
        # Reject vertices that are too far from the center
        max_reasonable_distance = min(enumerator.L ÷ 2, 4)  # Adaptive distance limit
        distance = lattice_distance(center_coords, new_coords, enumerator.L)
        
        return distance <= max_reasonable_distance
    end
    
    return true  # Allow all extensions for small loops
end

function is_valid_loop(vertices::Set{Int}, edges::Vector{Tuple{Int,Int}})
    if length(edges) < 2
        return false
    end
    
    degree = Dict{Int, Int}()
    for vertex in vertices
        degree[vertex] = 0
    end
    
    for (u, v) in edges
        degree[u] += 1
        degree[v] += 1
    end
    
    for vertex in vertices
        if degree[vertex] < 2
            return false
        end
    end
    
    return is_connected(vertices, edges)
end

function is_connected(vertices::Set{Int}, edges::Vector{Tuple{Int,Int}})
    if length(vertices) <= 1
        return true
    end
    
    adj = Dict{Int, Vector{Int}}()
    for vertex in vertices
        adj[vertex] = Int[]
    end
    
    for (u, v) in edges
        push!(adj[u], v)
        push!(adj[v], u)
    end
    
    start = first(vertices)
    visited = Set([start])
    queue = [start]
    
    while !isempty(queue)
        current = popfirst!(queue)
        for neighbor in adj[current]
            if neighbor in vertices && !(neighbor in visited)
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    
    return length(visited) == length(vertices)
end

function create_periodic_square_lattice(L::Int)
    N = L * L
    adj_matrix = zeros(Int, N, N)
    
    function coord_to_idx(i::Int, j::Int)
        return (i - 1) * L + j
    end
    
    for i in 1:L
        for j in 1:L
            current = coord_to_idx(i, j)
            
            right = coord_to_idx(i, j == L ? 1 : j + 1)
            adj_matrix[current, right] = 1
            adj_matrix[right, current] = 1
            
            down = coord_to_idx(i == L ? 1 : i + 1, j)
            adj_matrix[current, down] = 1
            adj_matrix[down, current] = 1
        end
    end
    
    return adj_matrix
end

function find_all_loops_in_graph(enumerator::LoopEnumerator, max_weight::Int)
    """Find all unique loops in the graph by enumerating from each vertex."""
    all_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    # Enumerate from each vertex and use canonical representation for deduplication
    for vertex in 1:enumerator.n_vertices
        vertex_loops = find_loops_supported_on_vertex(enumerator, vertex, max_weight)
        
        for loop in vertex_loops
            canonical = canonical_loop_representation(loop)
            
            if !(canonical in seen_loops)
                push!(seen_loops, canonical)
                new_loop = Loop(sort(loop.vertices), sort(loop.edges), loop.weight)
                push!(all_loops, new_loop)
            end
        end
    end
    
    return all_loops
end

function find_all_loops_in_graph_optimized(enumerator::LoopEnumerator, max_weight::Int)
    """Find all unique loops with speed optimizations while maintaining faithfulness."""
    all_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    # Optimization 1: Only enumerate from vertices that could be canonical representatives
    # This reduces redundant work while maintaining complete coverage
    for vertex in 1:enumerator.n_vertices
        vertex_loops = find_loops_supported_on_vertex_optimized(enumerator, vertex, max_weight, seen_loops)
        
        for loop in vertex_loops
            # Only keep loops where 'vertex' is the minimum vertex (canonical representative)
            if vertex == minimum(loop.vertices)
                canonical = canonical_loop_representation(loop)
                
                if !(canonical in seen_loops)
                    push!(seen_loops, canonical)
                    new_loop = Loop(sort(loop.vertices), sort(loop.edges), loop.weight)
                    push!(all_loops, new_loop)
                end
            end
        end
    end
    
    return all_loops
end