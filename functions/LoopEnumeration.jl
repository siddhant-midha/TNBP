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

struct LoopEnumerator
    adj_matrix::Matrix{Int}
    n_vertices::Int
    adj_list::Dict{Int, Vector{Int}}
    
    function LoopEnumerator(adj_matrix::Matrix{Int})
        n = size(adj_matrix, 1)
        adj_list = build_adjacency_list(adj_matrix, n)
        new(adj_matrix, n, adj_list)
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

function find_loops_supported_on_vertex(enumerator::LoopEnumerator, vertex::Int, max_weight::Int)
    found_loops = Loop[]
    find_loops_bfs!(enumerator, vertex, max_weight, found_loops)
    return found_loops
end

function find_loops_bfs!(enumerator::LoopEnumerator, start_vertex::Int, 
                        max_weight::Int, found_loops::Vector{Loop})
    neighbors = enumerator.adj_list[start_vertex]
    
    for i in 1:length(neighbors)
        for j in (i+1):length(neighbors)
            n1, n2 = neighbors[i], neighbors[j]
            initial_edges = [
                (min(start_vertex, n1), max(start_vertex, n1)),
                (min(start_vertex, n2), max(start_vertex, n2))
            ]
            initial_vertices = Set([start_vertex, n1, n2])
            
            if is_valid_loop(initial_vertices, initial_edges)
                push!(found_loops, Loop(sort(collect(initial_vertices)), sort(initial_edges), length(initial_edges)))
            end
            
            if length(initial_edges) < max_weight
                expand_subgraph!(enumerator, initial_vertices, initial_edges, 
                               max_weight, found_loops, start_vertex)
            end
        end
    end
end

function expand_subgraph!(enumerator::LoopEnumerator, vertices::Set{Int}, 
                         edges::Vector{Tuple{Int,Int}}, max_weight::Int,
                         found_loops::Vector{Loop}, target_vertex::Int)
    if length(edges) >= max_weight
        return
    end
    
    candidate_edges = Tuple{Int,Int}[]
    
    # Edges between existing vertices
    vertex_list = collect(vertices)
    for i in 1:length(vertex_list)
        for j in (i+1):length(vertex_list)
            v1, v2 = vertex_list[i], vertex_list[j]
            edge = (min(v1, v2), max(v1, v2))
            if !(edge in edges) && v2 in enumerator.adj_list[v1]
                push!(candidate_edges, edge)
            end
        end
    end
    
    # Edges to new vertices
    for v in vertices
        for neighbor in enumerator.adj_list[v]
            if !(neighbor in vertices)
                edge = (min(v, neighbor), max(v, neighbor))
                push!(candidate_edges, edge)
            end
        end
    end
    
    for edge in candidate_edges
        new_edges = copy(edges)
        push!(new_edges, edge)
        new_vertices = copy(vertices)
        union!(new_vertices, [edge[1], edge[2]])
        
        if !(target_vertex in new_vertices)
            continue
        end
        
        if is_valid_loop(new_vertices, new_edges)
            loop_candidate = Loop(sort(collect(new_vertices)), sort(new_edges), length(new_edges))
            
            is_duplicate = false
            for existing_loop in found_loops
                if (existing_loop.weight == loop_candidate.weight &&
                    existing_loop.vertices == loop_candidate.vertices &&
                    existing_loop.edges == loop_candidate.edges)
                    is_duplicate = true
                    break
                end
            end
            
            if !is_duplicate
                push!(found_loops, loop_candidate)
            end
        end
        
        if length(new_edges) < max_weight
            expand_subgraph!(enumerator, new_vertices, new_edges, max_weight, 
                           found_loops, target_vertex)
        end
    end
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

# No module exports needed for include-style usage