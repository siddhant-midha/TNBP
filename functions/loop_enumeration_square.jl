"""
loop_enumeration_square.jl

Efficient loop enumeration implementation for square lattices using coordinate-based
canonical forms during search and conversion to LoopEnumeration.jl format for output.
Supports both periodic and open boundary conditions.
"""

# Include the base LoopEnumeration.jl for data structures
include("LoopEnumeration.jl")

# Coordinate types from GPT implementation
const Coord = NTuple{2,Int}                       # (x,y) coordinates (0-based)
const Edge  = NTuple{2,Coord}                     # ((x1,y1),(x2,y2))

# Enhanced LoopEnumerator with boundary condition support
struct SquareLoopEnumerator
    base_enumerator::LoopEnumerator
    L::Int                          # Lattice size
    periodic::Bool                  # Boundary condition flag
    coord_to_vertex::Dict{Coord, Int}  # (x,y) -> vertex index
    vertex_to_coord::Dict{Int, Coord}  # vertex index -> (x,y)
    
    function SquareLoopEnumerator(L::Int; periodic::Bool=true)
        adj_matrix = periodic ? create_periodic_square_lattice(L) : create_open_square_lattice(L)
        base_enum = LoopEnumerator(adj_matrix)
        
        # Build coordinate mappings
        coord_to_vertex = Dict{Coord, Int}()
        vertex_to_coord = Dict{Int, Coord}()
        
        for i in 0:L-1, j in 0:L-1
            coord = (i, j)
            vertex = coord_to_vertex_index(i, j, L)  # 1-based vertex index
            coord_to_vertex[coord] = vertex
            vertex_to_coord[vertex] = coord
        end
        
        new(base_enum, L, periodic, coord_to_vertex, vertex_to_coord)
    end
end

# ========== Coordinate Conversion Utilities ==========

"""Convert (i,j) coordinates to vertex index (1-based Julia indexing)"""
@inline function coord_to_vertex_index(i::Int, j::Int, L::Int)
    return i * L + j + 1
end

"""Convert vertex index to (i,j) coordinates (0-based)"""
@inline function vertex_index_to_coord(vertex::Int, L::Int)
    vertex_0based = vertex - 1
    i = vertex_0based ÷ L
    j = vertex_0based % L
    return (i, j)
end

"""Create open boundary square lattice adjacency matrix"""
function create_open_square_lattice(L::Int)
    N = L * L
    adj_matrix = zeros(Int, N, N)
    
    for i in 0:L-1, j in 0:L-1
        current = coord_to_vertex_index(i, j, L)
        
        # Right neighbor (within bounds only)
        if j < L-1
            right = coord_to_vertex_index(i, j+1, L)
            adj_matrix[current, right] = 1
            adj_matrix[right, current] = 1
        end
        
        # Down neighbor (within bounds only)  
        if i < L-1
            down = coord_to_vertex_index(i+1, j, L)
            adj_matrix[current, down] = 1
            adj_matrix[down, current] = 1
        end
    end
    
    return adj_matrix
end
# ========== D4 Symmetry Transformations ==========

"""Apply coordinate transformation with boundary condition awareness"""
@inline function apply_coord_transform(coord::Coord, transform_func::Function, L::Int, periodic::Bool)
    if periodic
        # Apply transformation and wrap coordinates
        new_coord = transform_func(coord)
        x_wrapped = mod(new_coord[1], L)
        y_wrapped = mod(new_coord[2], L)
        return (x_wrapped, y_wrapped)
    else
        # Apply transformation and check bounds
        new_coord = transform_func(coord)
        if 0 <= new_coord[1] < L && 0 <= new_coord[2] < L
            return new_coord
        else
            return nothing  # Out of bounds for open boundary
        end
    end
end

# D4 group transformations (same as GPT version but with boundary awareness)
@inline T_id(p::Coord)      = ( p[1],  p[2])
@inline T_r90(p::Coord)     = (-p[2],  p[1])
@inline T_r180(p::Coord)    = (-p[1], -p[2])
@inline T_r270(p::Coord)    = ( p[2], -p[1])
@inline T_refx(p::Coord)    = (-p[1],  p[2])    # reflect across y-axis
@inline T_refy(p::Coord)    = ( p[1], -p[2])    # reflect across x-axis
@inline T_refdiag(p::Coord) = ( p[2],  p[1])    # reflect across y=x
@inline T_refanti(p::Coord) = (-p[2], -p[1])    # reflect across y=-x

const D4 = (T_id, T_r90, T_r180, T_r270, T_refx, T_refy, T_refdiag, T_refanti)

# ========== Edge and Canonical Form Utilities ==========

"""Lexicographic order for coordinates"""
@inline leq_coord(a::Coord, b::Coord) = (a[1] < b[1]) || (a[1] == b[1] && a[2] <= b[2])

"""Normalized, ordered edge"""
@inline norm_edge(u::Coord, v::Coord)::Edge = leq_coord(u,v) ? (u,v) : (v,u)

"""Convert edge set to sortable string key"""
function edge_list_key(E::Set{Edge})
    es = collect(E)
    sort!(es, by = e -> (e[1][1], e[1][2], e[2][1], e[2][2]))
    io = IOBuffer()
    for e in es
        print(io, e[1][1], ',', e[1][2], ';', e[2][1], ',', e[2][2], '|')
    end
    return String(take!(io))
end

"""Canonical key under D4 symmetries (GPT-style for search efficiency)"""
function gpt_canonical_key(E::Set{Edge}, L::Int, periodic::Bool)
    best::Union{Nothing,String} = nothing
    for T in D4
        # Transform and renormalize each edge
        E2 = Set{Edge}()
        valid = true
        for (u,v) in E
            u_new = apply_coord_transform(u, T, L, periodic)
            v_new = apply_coord_transform(v, T, L, periodic)
            
            if u_new === nothing || v_new === nothing
                valid = false
                break
            end
            
            push!(E2, norm_edge(u_new, v_new))
        end
        
        if valid
            k = edge_list_key(E2)
            if best === nothing || k < best
                best = k
            end
        end
    end
    return best
end

"""Get canonical edge set under D4 (for transformation recovery)"""
function gpt_canonical_edges(E::Set{Edge}, L::Int, periodic::Bool)
    best_key::Union{Nothing,String} = nothing
    best_E::Union{Nothing,Set{Edge}} = nothing
    
    for T in D4
        E2 = Set{Edge}()
        valid = true
        for (u,v) in E
            u_new = apply_coord_transform(u, T, L, periodic)
            v_new = apply_coord_transform(v, T, L, periodic)
            
            if u_new === nothing || v_new === nothing
                valid = false
                break
            end
            
            push!(E2, norm_edge(u_new, v_new))
        end
        
        if valid
            k = edge_list_key(E2)
            if best_key === nothing || k < best_key
                best_key, best_E = k, E2
            end
        end
    end
    return best_E
end

# ========== BFS Loop Enumeration State ==========

mutable struct SquareLoopState
    V::Set{Coord}               # vertices present (in coordinate form)
    E::Set{Edge}                # edges present (in coordinate form)  
    deg::Dict{Coord,Int}        # current degrees within the subgraph
end

"""Get neighbors by directly checking adjacency matrix"""
function get_coord_neighbors(v::Coord, L::Int, periodic::Bool)
    neighbors = Coord[]
    v_idx = coord_to_vertex_index(v[1], v[2], L)
    
    # Get adjacency matrix
    adj_matrix = periodic ? create_periodic_square_lattice(L) : create_open_square_lattice(L)
    
    # Check all possible vertices
    for u_idx in 1:(L*L)
        if adj_matrix[v_idx, u_idx] == 1
            u_coord = vertex_index_to_coord(u_idx, L)
            push!(neighbors, u_coord)
        end
    end
    
    return neighbors
end


"""Compute frontier edges: all lattice edges with ≥1 endpoint in V that are not yet in E"""
function frontier_edges(V::Set{Coord}, E::Set{Edge}, L::Int, periodic::Bool)
    F = Set{Edge}()
    for v in V
        for u in get_coord_neighbors(v, L, periodic)
            e = norm_edge(v, u)
            if !(e in E)
                # For open boundaries, double-check that this is a valid lattice edge
                if !periodic
                    # Only add edges between immediate grid neighbors (not wrapping)
                    dx = abs(v[1] - u[1])
                    dy = abs(v[2] - u[2])
                    # Valid edge: Manhattan distance = 1 and no coordinate wrapping
                    if (dx == 1 && dy == 0) || (dx == 0 && dy == 1)
                        push!(F, e)
                    end
                else
                    push!(F, e)
                end
            end
        end
    end
    return F
end
"""Add an edge and update degrees/vertex set"""
function add_edge(s::SquareLoopState, e::Edge)
    (u,v) = e
    V2 = copy(s.V)
    E2 = copy(s.E)
    deg2 = copy(s.deg)
    push!(E2, e)
    if !(u in V2); push!(V2,u); deg2[u] = 0; end
    if !(v in V2); push!(V2,v); deg2[v] = 0; end
    deg2[u] = get(deg2,u,0) + 1
    deg2[v] = get(deg2,v,0) + 1
    return SquareLoopState(V2, E2, deg2)
end

"""Sum of degree deficits"""
@inline total_deficit(deg::Dict{Coord,Int}) = sum(d < 2 ? 2 - d : 0 for d in values(deg))

# ========== Main Loop Enumeration Function ==========

"""
Enumerate loops supported on a given coordinate using GPT-style BFS with D4 canonical forms
"""
function enumerate_loops_on_coord(root_coord::Coord, max_edges::Int, L::Int, periodic::Bool; min_edges::Int=2)
    # Initial state: only the root present, no edges, degree 0
    init = SquareLoopState(Set([root_coord]), Set{Edge}(), Dict(root_coord => 0))

    # BFS queue and head pointer
    queue = Vector{SquareLoopState}([init])
    head = 1

    # Seen canonical keys per edge-count
    seen = Dict{Int, Set{String}}()
    seen[0] = Set([gpt_canonical_key(init.E, L, periodic)])

    # Collected results by size
    results = Dict{Int, Vector{SquareLoopState}}()

    while head <= length(queue)
        s = queue[head]; head += 1
        m = length(s.E)

        # Get canonical key for duplicate detection
        canonical_key = gpt_canonical_key(s.E, L, periodic)
        if canonical_key === nothing
            continue
        end
        

        # Early pruning
        deficit = total_deficit(s.deg)
        remaining = max_edges - m
        if remaining < cld(deficit, 2)   # ceil(deficit/2)
            continue
        end

        # If this is a valid loop, record it
        if deficit == 0 && m >= min_edges
            push!(get!(results, m, Vector{SquareLoopState}()), s)
        end

        # Do not expand beyond max_edges
        if m == max_edges
            continue
        end

        # Expand by adding each frontier edge
        for e in frontier_edges(s.V, s.E, L, periodic)
            s2 = add_edge(s, e)
            m2 = m + 1
            cank = gpt_canonical_key(s2.E, L, periodic)
            if cank !== nothing
                bucket = get!(seen, m2, Set{String}())
                if !(cank in bucket)
                    push!(bucket, cank)
                    push!(queue, s2)
                end
            end
        end
    end

    return results
end

"""
Convert coordinate-based loop states to Loop objects and apply all D4 transformations
"""
function convert_and_expand_loops(coord_results::Dict{Int, Vector{SquareLoopState}}, 
                                  enumerator::SquareLoopEnumerator)
    all_loops = Loop[]
    seen_canonical = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    for (weight, states) in coord_results
        for state in states
            # Get the canonical edge set under D4
            canonical_edges = gpt_canonical_edges(state.E, enumerator.L, enumerator.periodic)
            if canonical_edges === nothing
                continue
            end
            
            # Apply all D4 transformations to recover equivalent loops
            for T in D4
                transformed_edges = Set{Edge}()
                valid = true
                
                for (u, v) in canonical_edges
                    u_new = apply_coord_transform(u, T, enumerator.L, enumerator.periodic)
                    v_new = apply_coord_transform(v, T, enumerator.L, enumerator.periodic)
                    
                    if u_new === nothing || v_new === nothing
                        valid = false
                        break
                    end
                    
                    push!(transformed_edges, norm_edge(u_new, v_new))
                end
                
                if !valid
                    continue
                end
                
                # Convert back to vertex indices
                vertices = Set{Int}()
                edges_vertex = Tuple{Int,Int}[]
                
                for (u_coord, v_coord) in transformed_edges
                    u_vertex = enumerator.coord_to_vertex[u_coord]
                    v_vertex = enumerator.coord_to_vertex[v_coord]
                    push!(vertices, u_vertex)
                    push!(vertices, v_vertex)
                    push!(edges_vertex, (min(u_vertex, v_vertex), max(u_vertex, v_vertex)))
                end
                
                # Create Loop object
                loop = Loop(sort(collect(vertices)), sort(edges_vertex), weight)
                canonical = canonical_loop_representation(loop)
                
                if !(canonical in seen_canonical)
                    push!(seen_canonical, canonical)
                    push!(all_loops, loop)
                end
            end
        end
    end
    
    return all_loops
end

# ========== Main Interface Functions ==========

"""
Find all loops supported on a vertex using enhanced square lattice enumeration
"""
function find_loops_supported_on_vertex_square(enumerator::SquareLoopEnumerator, vertex::Int, max_weight::Int)
    # Convert vertex to coordinate
    root_coord = enumerator.vertex_to_coord[vertex]
    
    # Enumerate loops in coordinate space with D4 canonical forms
    coord_results = enumerate_loops_on_coord(root_coord, max_weight, enumerator.L, enumerator.periodic)
    
    # Convert back and expand with all transformations
    loops = convert_and_expand_loops(coord_results, enumerator)
    
    # Filter to only include loops that actually contain the target vertex
    filtered_loops = Loop[]
    for loop in loops
        if vertex in loop.vertices
            push!(filtered_loops, loop)
        end
    end
    
    return filtered_loops
end

# Compatibility wrapper removed - use original LoopEnumeration.jl function instead

"""
Find all unique loops in the graph using square lattice enumeration
"""
function find_all_loops_in_graph_square(enumerator::SquareLoopEnumerator, max_weight::Int)
    all_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    # Only enumerate from canonical representative vertices to avoid duplicates
    for vertex in 1:enumerator.base_enumerator.n_vertices
        vertex_loops = find_loops_supported_on_vertex_square(enumerator, vertex, max_weight)
        
        for loop in vertex_loops
            # Only keep loops where 'vertex' is the minimum vertex
            if vertex == minimum(loop.vertices)
                canonical = canonical_loop_representation(loop)
                
                if !(canonical in seen_loops)
                    push!(seen_loops, canonical)
                    push!(all_loops, loop)
                end
            end
        end
    end
    
    return all_loops
end

# Compatibility wrapper removed - use original LoopEnumeration.jl function instead