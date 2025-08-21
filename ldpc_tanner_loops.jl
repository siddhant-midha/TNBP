"""
LDPC Tanner Graph Loop Enumeration

This module provides functions to find loops in the Tanner graph of an LDPC code
that are supported on a given data bit.

Author: Generated for TNBP project
"""

using SparseArrays

"""
    TannerLoop

Represents a loop in the Tanner graph of an LDPC code.

Fields:
- `edges`: Vector of edges in the form [(v1,v2), (v2,v3), ..., (vn,v1)]
- `data_bits`: Set of data bit indices involved in the loop
- `check_bits`: Set of check bit indices involved in the loop
- `length`: Number of edges in the loop
"""
struct TannerLoop
    edges::Vector{Tuple{Int,Int}}
    data_bits::Set{Int}
    check_bits::Set{Int}
    length::Int
end

"""
    get_tanner_neighbors(H::AbstractMatrix, bit_index::Int, is_data_bit::Bool)

Get neighbors of a bit in the Tanner graph.

Arguments:
- `H`: Parity check matrix (m × n)
- `bit_index`: Index of the bit
- `is_data_bit`: true if it's a data bit (column), false if check bit (row)

Returns:
- Vector of neighbor indices
"""
function get_tanner_neighbors(H::AbstractMatrix, bit_index::Int, is_data_bit::Bool)
    if is_data_bit
        # Data bit: find which check bits it's connected to
        return findall(x -> x != 0, H[:, bit_index])
    else
        # Check bit: find which data bits it's connected to
        return findall(x -> x != 0, H[bit_index, :])
    end
end

"""
    bipartite_bfs_loops(H::AbstractMatrix, start_data_bit::Int, max_length::Int=20)

Find all loops in the Tanner graph supported on a given data bit using BFS.

Arguments:
- `H`: Parity check matrix (m × n)
- `start_data_bit`: Data bit index (1 to n)
- `max_length`: Maximum loop length to search

Returns:
- Vector of TannerLoop objects
"""
function bipartite_bfs_loops(H::AbstractMatrix, start_data_bit::Int, max_length::Int=20)
    m, n = size(H)
    
    # Validate input
    if start_data_bit < 1 || start_data_bit > n
        throw(ArgumentError("start_data_bit must be between 1 and $n"))
    end
    
    loops = TannerLoop[]
    
    # Get initial check bit neighbors of the starting data bit
    initial_checks = get_tanner_neighbors(H, start_data_bit, true)
    
    # For each pair of initial check bits, try to find paths back to start_data_bit
    for i in 1:length(initial_checks)
        for j in (i+1):length(initial_checks)
            check1, check2 = initial_checks[i], initial_checks[j]
            
            # Find all paths from check1 to check2 that could form loops
            paths = find_bipartite_paths(H, check1, check2, start_data_bit, max_length÷2)
            
            for path in paths
                # Construct the complete loop
                loop_edges = construct_loop_edges(start_data_bit, check1, check2, path, n)
                if !isempty(loop_edges)
                    data_bits, check_bits = extract_bits_from_edges(loop_edges, n)
                    
                    # Only include if the starting data bit is in the loop
                    if start_data_bit in data_bits
                        loop = TannerLoop(loop_edges, data_bits, check_bits, length(loop_edges))
                        push!(loops, loop)
                    end
                end
            end
        end
    end
    
    # Remove duplicate loops
    unique_loops = remove_duplicate_loops(loops)
    
    return unique_loops
end

"""
    find_bipartite_paths(H::AbstractMatrix, start_check::Int, end_check::Int, 
                        avoid_data::Int, max_steps::Int)

Find paths between two check bits in the bipartite Tanner graph.
"""
function find_bipartite_paths(H::AbstractMatrix, start_check::Int, end_check::Int, 
                             avoid_data::Int, max_steps::Int)
    m, n = size(H)
    paths = Vector{Vector{Int}}()
    
    # BFS to find paths
    # Each state is (current_node, is_data_bit, path, steps)
    queue = [(start_check, false, [start_check], 0)]
    
    while !isempty(queue)
        current, is_data, path, steps = popfirst!(queue)
        
        if steps >= max_steps
            continue
        end
        
        if is_data
            # Currently at a data bit, go to check bits
            neighbors = get_tanner_neighbors(H, current, true)
            for neighbor in neighbors
                if neighbor == end_check && steps > 0
                    # Found a path to the end check
                    push!(paths, vcat(path, [neighbor]))
                elseif neighbor != end_check && neighbor != start_check
                    # Continue exploring
                    new_path = vcat(path, [neighbor])
                    push!(queue, (neighbor, false, new_path, steps + 1))
                end
            end
        else
            # Currently at a check bit, go to data bits
            neighbors = get_tanner_neighbors(H, current, false)
            for neighbor in neighbors
                if neighbor != avoid_data && !(neighbor in path)
                    new_path = vcat(path, [neighbor])
                    push!(queue, (neighbor, true, new_path, steps + 1))
                end
            end
        end
    end
    
    return paths
end

"""
    ordered_edge(v1::Int, v2::Int)

Create an edge tuple with the smaller vertex first, ensuring v1 < v2.
"""
function ordered_edge(v1::Int, v2::Int)
    return v1 < v2 ? (v1, v2) : (v2, v1)
end

"""
    construct_loop_edges(start_data::Int, check1::Int, check2::Int, 
                        path::Vector{Int}, n::Int)

Construct the edge list for a loop in the Tanner graph.
"""
function construct_loop_edges(start_data::Int, check1::Int, check2::Int, 
                             path::Vector{Int}, n::Int)
    edges = Tuple{Int,Int}[]
    
    # Add edge from start_data to check1
    push!(edges, ordered_edge(start_data, n + check1))
    
    # Add edges along the path
    for idx in 1:(length(path)-1)
        v1, v2 = path[idx], path[idx+1]
        # Determine if nodes are data or check bits based on the alternating pattern
        if idx % 2 == 1  # v1 is check, v2 is data
            push!(edges, ordered_edge(n + v1, v2))
        else  # v1 is data, v2 is check
            push!(edges, ordered_edge(v1, n + v2))
        end
    end
    
    # Add edge from check2 back to start_data
    push!(edges, ordered_edge(n + check2, start_data))
    
    return edges
end

"""
    extract_bits_from_edges(edges::Vector{Tuple{Int,Int}}, n::Int)

Extract data and check bit sets from edge list.
"""
function extract_bits_from_edges(edges::Vector{Tuple{Int,Int}}, n::Int)
    data_bits = Set{Int}()
    check_bits = Set{Int}()
    
    for (v1, v2) in edges
        if v1 <= n
            push!(data_bits, v1)
        else
            push!(check_bits, v1 - n)
        end
        
        if v2 <= n
            push!(data_bits, v2)
        else
            push!(check_bits, v2 - n)
        end
    end
    
    return data_bits, check_bits
end

"""
    remove_duplicate_loops(loops::Vector{TannerLoop})

Remove duplicate loops from the list.
"""
function remove_duplicate_loops(loops::Vector{TannerLoop})
    unique_loops = TannerLoop[]
    seen_signatures = Set{Tuple{Set{Int}, Set{Int}}}()
    
    for loop in loops
        signature = (loop.data_bits, loop.check_bits)
        if !(signature in seen_signatures)
            push!(seen_signatures, signature)
            push!(unique_loops, loop)
        end
    end
    
    return unique_loops
end

"""
    find_all_cycles_dfs(H::AbstractMatrix, start_data_bit::Int, max_length::Int=20)

Find all cycles in the Tanner graph using DFS approach.
This is an alternative implementation using depth-first search.
"""
function find_all_cycles_dfs(H::AbstractMatrix, start_data_bit::Int, max_length::Int=20)
    m, n = size(H)
    
    if start_data_bit < 1 || start_data_bit > n
        throw(ArgumentError("start_data_bit must be between 1 and $n"))
    end
    
    cycles = TannerLoop[]
    visited = Set{Tuple{Int,Bool}}  # (node_index, is_data_bit)
    
    # Start DFS from the given data bit
    dfs_find_cycles!(H, start_data_bit, true, start_data_bit, 
                     Tuple{Int,Int}[], visited, cycles, max_length, n)
    
    return remove_duplicate_loops(cycles)
end

"""
    dfs_find_cycles!(H, current, is_data, start, path, visited, cycles, max_len, n)

Recursive DFS helper for finding cycles.
"""
function dfs_find_cycles!(H::AbstractMatrix, current::Int, is_data::Bool, 
                         start::Int, path::Vector{Tuple{Int,Int}}, 
                         visited::Set{Tuple{Int,Bool}}, cycles::Vector{TannerLoop},
                         max_len::Int, n::Int)
    
    if length(path) > max_len
        return
    end
    
    # Get neighbors
    neighbors = get_tanner_neighbors(H, current, is_data)
    
    for neighbor in neighbors
        next_is_data = !is_data
        state = (neighbor, next_is_data)
        
        # Check if we've found a cycle back to start
        if neighbor == start && next_is_data && length(path) > 2
            # Complete the cycle
            final_edge = is_data ? ordered_edge(current, n + neighbor) : ordered_edge(n + current, neighbor)
            cycle_edges = vcat(path, [final_edge])
            
            data_bits, check_bits = extract_bits_from_edges(cycle_edges, n)
            cycle = TannerLoop(cycle_edges, data_bits, check_bits, length(cycle_edges))
            push!(cycles, cycle)
            continue
        end
        
        # Continue DFS if we haven't visited this state and it's not the start
        if !(state in visited) && !(neighbor == start && length(path) <= 2)
            push!(visited, state)
            
            new_edge = is_data ? ordered_edge(current, n + neighbor) : ordered_edge(n + current, neighbor)
            new_path = vcat(path, [new_edge])
            
            dfs_find_cycles!(H, neighbor, next_is_data, start, new_path, 
                           visited, cycles, max_len, n)
            
            delete!(visited, state)
        end
    end
end

"""
    print_loop(loop::TannerLoop)

Pretty print a loop for debugging.
"""
function print_loop(loop::TannerLoop)
    println("Loop with $(loop.length) edges:")
    println("  Edges: ", loop.edges)
    println("  Data bits: ", sort(collect(loop.data_bits)))
    println("  Check bits: ", sort(collect(loop.check_bits)))
end

"""
    find_tanner_loops(H::AbstractMatrix, data_bit_index::Int; 
                      max_length::Int=20, method::Symbol=:bfs)

Main function to find loops in the Tanner graph supported on a given data bit.

Arguments:
- `H`: Parity check matrix (m × n)
- `data_bit_index`: Data bit index (1 to n) 
- `max_length`: Maximum loop length to search (default: 20)
- `method`: Search method (:efficient, :bfs, or :dfs, default: :efficient)

Returns:
- Vector of TannerLoop objects representing loops supported on the data bit
- Each loop contains edges in the format [(v1,v2), (v2,v3), ..., (vn,v1)]
- All edges follow the convention v1 < v2 (smaller vertex index first)
- Data bits are numbered 1:n, check bits are numbered (n+1):(n+m)

The :efficient method finds all loops in the graph once and filters for the given data bit.
This ensures that all loops containing the data bit are found, not just those starting from it.

Example:
```julia
# Create a simple LDPC parity check matrix
H = [1 1 0 1 0 0;
     0 1 1 0 1 0;
     1 0 1 0 0 1]

# Find loops supported on data bit 1
loops = find_tanner_loops(H, 1)

# Print all loops
for (i, loop) in enumerate(loops)
    println("Loop ", i, ":")
    print_loop(loop)
    println()
end
```
"""
function find_tanner_loops(H::AbstractMatrix, data_bit_index::Int; 
                          max_length::Int=20, method::Symbol=:efficient)
    
    if method == :efficient
        # Find all loops in the graph and filter for ones containing the data bit
        all_loops = find_all_loops_in_tanner_graph(H, max_length)
        return filter(loop -> data_bit_index in loop.data_bits, all_loops)
    elseif method == :bfs
        return bipartite_bfs_loops(H, data_bit_index, max_length)
    elseif method == :dfs
        return find_all_cycles_dfs(H, data_bit_index, max_length)
    else
        throw(ArgumentError("method must be :efficient, :bfs, or :dfs"))
    end
end

"""
    find_all_loops_in_tanner_graph(H::AbstractMatrix, max_length::Int=20)

Find ALL loops in the Tanner graph, regardless of which data bit they start from.

Arguments:
- `H`: Parity check matrix (m × n)
- `max_length`: Maximum loop length to search (default: 20)

Returns:
- Vector of TannerLoop objects representing all loops in the graph
"""
function find_all_loops_in_tanner_graph(H::AbstractMatrix, max_length::Int=20)
    m, n = size(H)
    all_loops = TannerLoop[]
    
    # Find loops starting from each data bit and collect them all
    for data_bit in 1:n
        loops = bipartite_bfs_loops(H, data_bit, max_length)
        append!(all_loops, loops)
    end
    
    # Remove duplicates based on the set of vertices involved
    return remove_duplicate_loops(all_loops)
end

"""
    validate_loop(edges)

Validate that a loop is actually a valid cycle.
"""
function validate_loop(edges)
    if length(edges) < 3
        return false
    end
    
    # Check that edges form a proper cycle
    for i in 1:length(edges)
        current_edge = edges[i]
        next_edge = edges[i % length(edges) + 1]
        
        # The second vertex of current edge should be the first vertex of next edge
        if current_edge[2] != next_edge[1]
            return false
        end
    end
    
    # Check no repeated edges
    if length(unique(edges)) != length(edges)
        return false
    end
    
    return true
end

"""
    create_simple_example()

Create a simple LDPC matrix example for demonstration.
"""
function create_simple_example()
    # Simple 3x6 parity check matrix
    H = [1 1 0 1 0 0;
         0 1 1 0 1 0;
         1 0 1 0 0 1]
    
    println("Simple LDPC Example")
    println("===================")
    println("Parity check matrix H:")
    display(H)
    println()
    
    # Find all loops for data bit 1
    loops = find_tanner_loops(H, 1)
    println("Loops supported on data bit 1:")
    
    if isempty(loops)
        println("  No loops found")
    else
        for (i, loop) in enumerate(loops)
            println("  Loop $i (length $(loop.length)):")
            println("    Edges: ", loop.edges)
            println("    Data bits involved: ", sort(collect(loop.data_bits)))
            println("    Check bits involved: ", sort(collect(loop.check_bits)))
            println()
        end
    end
    
    return H, loops
end

# Example usage and test function
"""
    test_tanner_loops()

Test the loop finding functionality with a simple example.
"""
function test_tanner_loops()
    println("Testing Tanner Graph Loop Enumeration")
    println("="^50)
    
    # Create a simple LDPC parity check matrix
    H = [1 1 0 1 0 0;
         0 1 1 0 1 0;
         1 0 1 0 0 1]
    
    println("Parity check matrix H:")
    display(H)
    println()
    
    m, n = size(H)
    println("Matrix dimensions: $m × $n (m=$m check bits, n=$n data bits)")
    println()
    
    # Test for each data bit
    for data_bit in 1:n
        println("Finding loops supported on data bit $data_bit:")
        loops = find_tanner_loops(H, data_bit, max_length=10)
        
        if isempty(loops)
            println("  No loops found.")
        else
            for (i, loop) in enumerate(loops)
                println("  Loop $i:")
                println("    Edges: ", loop.edges)
                println("    Data bits: ", sort(collect(loop.data_bits)))
                println("    Check bits: ", sort(collect(loop.check_bits)))
            end
        end
        println()
    end
end

"""
    test_tanner_loops_comprehensive()

More comprehensive test with a matrix that should have multiple loops.
"""
function test_tanner_loops_comprehensive()
    println("Comprehensive Test: Tanner Graph Loop Enumeration")
    println("=======================================================")
    
    # Test matrix with known 4-cycles
    H = [1 1 1 0 0 0;
         1 0 0 1 1 0;
         0 1 0 1 0 1;
         0 0 1 0 1 1]
    
    println("Parity check matrix H:")
    display(H)
    println("Matrix dimensions: $(size(H, 1)) × $(size(H, 2)) (m=$(size(H, 1)) check bits, n=$(size(H, 2)) data bits)")
    println()
    
    # Test all data bits
    for data_bit in 1:size(H, 2)
        println("Finding loops supported on data bit $data_bit:")
        loops = find_tanner_loops(H, data_bit)
        
        if isempty(loops)
            println("  No loops found.")
        else
            for (i, loop) in enumerate(loops)
                println("  Loop $i:")
                println("    Edges: ", loop.edges)
                println("    Data bits: ", sort(collect(loop.data_bits)))
                println("    Check bits: ", sort(collect(loop.check_bits)))
            end
        end
        println()
    end
    
    # Also test finding shortest cycle containing a specific data bit
    println("Testing shortest loops:")
    for data_bit in 1:size(H, 2)
        loops = find_tanner_loops(H, data_bit, max_length=8)
        if !isempty(loops)
            shortest = minimum(loop.length for loop in loops)
            println("Data bit $data_bit: shortest loop length = $shortest")
        else
            println("Data bit $data_bit: no loops found")
        end
    end
end

"""
    get_loops_for_data_bit(H, data_bit_index; max_length=20)

Main API function to get loops in the Tanner graph supported on a given data bit.

Arguments:
- `H`: Parity check matrix (m × n) 
- `data_bit_index`: Data bit index (1 to n)
- `max_length`: Maximum loop length to search (default: 20)

Returns:
- Vector of vectors, where each inner vector represents a loop as a list of edges
  Each edge is a tuple (v1, v2) where:
  - Data bits are numbered 1:n  
  - Check bits are numbered (n+1):(n+m)
  - Loop is represented as [(v1,v2), (v2,v3), ..., (vn,v1)]

Example:
```julia
H = [1 1 0; 1 0 1; 0 1 1]  # 3x3 parity check matrix
loops = get_loops_for_data_bit(H, 1)
# Returns: [[(1, 4), (4, 2), (2, 5), (5, 3), (3, 6), (6, 1)]]
```
"""
function get_loops_for_data_bit(H::AbstractMatrix, data_bit_index::Int; max_length::Int=20)
    # Find all TannerLoop objects
    tanner_loops = find_tanner_loops(H, data_bit_index, max_length=max_length)
    
    # Convert to simple edge lists
    return [loop.edges for loop in tanner_loops]
end

"""
    find_loops_for_check_bit(H::AbstractMatrix, check_bit_index::Int; 
                            max_length::Int=20, method::Symbol=:efficient)

Find loops in the Tanner graph supported on a given check bit.

Arguments:
- `H`: Parity check matrix (m × n)
- `check_bit_index`: Check bit index (1 to m) 
- `max_length`: Maximum loop length to search (default: 20)
- `method`: Search method (:efficient or :bfs, default: :efficient)

Returns:
- Vector of TannerLoop objects representing loops supported on the check bit
- Each loop contains edges in the format [(v1,v2), (v2,v3), ..., (vn,v1)]
- All edges follow the convention v1 < v2 (smaller vertex index first)
- Data bits are numbered 1:n, check bits are numbered (n+1):(n+m)

Example:
```julia
# Create a simple LDPC parity check matrix
H = [1 1 0 1 0 0;
     0 1 1 0 1 0;
     1 0 1 0 0 1]

# Find loops supported on check bit 1
loops = find_loops_for_check_bit(H, 1)

# Print all loops
for (i, loop) in enumerate(loops)
    println("Loop ", i, ":")
    print_loop(loop)
    println()
end
```
"""
function find_loops_for_check_bit(H::AbstractMatrix, check_bit_index::Int; 
                                 max_length::Int=20, method::Symbol=:efficient)
    m, n = size(H)
    
    if check_bit_index < 1 || check_bit_index > m
        throw(ArgumentError("check_bit_index must be between 1 and $m"))
    end
    
    if method == :efficient
        # Find all loops in the graph and filter for ones containing the check bit
        all_loops = find_all_loops_in_tanner_graph(H, max_length)
        return filter(loop -> check_bit_index in loop.check_bits, all_loops)
    elseif method == :bfs
        return bipartite_bfs_loops_check(H, check_bit_index, max_length)
    else
        throw(ArgumentError("method must be :efficient or :bfs for check bit searches"))
    end
end

"""
    bipartite_bfs_loops_check(H::AbstractMatrix, start_check_bit::Int, max_length::Int=20)

Find all loops in the Tanner graph supported on a given check bit using BFS.

Arguments:
- `H`: Parity check matrix (m × n)
- `start_check_bit`: Check bit index (1 to m)
- `max_length`: Maximum loop length to search

Returns:
- Vector of TannerLoop objects
"""
function bipartite_bfs_loops_check(H::AbstractMatrix, start_check_bit::Int, max_length::Int=20)
    m, n = size(H)
    
    # Validate input
    if start_check_bit < 1 || start_check_bit > m
        throw(ArgumentError("start_check_bit must be between 1 and $m"))
    end
    
    loops = TannerLoop[]
    
    # Get initial data bit neighbors of the starting check bit
    initial_data = get_tanner_neighbors(H, start_check_bit, false)
    
    # For each pair of initial data bits, try to find paths back to start_check_bit
    for i in 1:length(initial_data)
        for j in (i+1):length(initial_data)
            data1, data2 = initial_data[i], initial_data[j]
            
            # Find all paths from data1 to data2 that could form loops
            paths = find_bipartite_paths_check(H, data1, data2, start_check_bit, max_length÷2)
            
            for path in paths
                # Construct the complete loop
                loop_edges = construct_loop_edges_check(start_check_bit, data1, data2, path, n)
                if !isempty(loop_edges)
                    data_bits, check_bits = extract_bits_from_edges(loop_edges, n)
                    
                    # Only include if the starting check bit is in the loop
                    if start_check_bit in check_bits
                        loop = TannerLoop(loop_edges, data_bits, check_bits, length(loop_edges))
                        push!(loops, loop)
                    end
                end
            end
        end
    end
    
    # Remove duplicate loops
    unique_loops = remove_duplicate_loops(loops)
    
    return unique_loops
end

"""
    find_bipartite_paths_check(H::AbstractMatrix, start_data::Int, end_data::Int, 
                              avoid_check::Int, max_steps::Int)

Find paths between two data bits in the bipartite Tanner graph.
"""
function find_bipartite_paths_check(H::AbstractMatrix, start_data::Int, end_data::Int, 
                                   avoid_check::Int, max_steps::Int)
    m, n = size(H)
    paths = Vector{Vector{Int}}()
    
    # BFS to find paths
    # Each state is (current_node, is_data_bit, path, steps)
    queue = [(start_data, true, [start_data], 0)]
    
    while !isempty(queue)
        current, is_data, path, steps = popfirst!(queue)
        
        if steps >= max_steps
            continue
        end
        
        if is_data
            # Currently at a data bit, go to check bits
            neighbors = get_tanner_neighbors(H, current, true)
            for neighbor in neighbors
                if neighbor != avoid_check && !(neighbor in path)
                    new_path = vcat(path, [neighbor])
                    push!(queue, (neighbor, false, new_path, steps + 1))
                end
            end
        else
            # Currently at a check bit, go to data bits
            neighbors = get_tanner_neighbors(H, current, false)
            for neighbor in neighbors
                if neighbor == end_data && steps > 0
                    # Found a path to the end data bit
                    push!(paths, vcat(path, [neighbor]))
                elseif neighbor != end_data && neighbor != start_data
                    # Continue exploring
                    new_path = vcat(path, [neighbor])
                    push!(queue, (neighbor, true, new_path, steps + 1))
                end
            end
        end
    end
    
    return paths
end

"""
    construct_loop_edges_check(start_check::Int, data1::Int, data2::Int, 
                              path::Vector{Int}, n::Int)

Construct the edge list for a loop starting from a check bit in the Tanner graph.
"""
function construct_loop_edges_check(start_check::Int, data1::Int, data2::Int, 
                                   path::Vector{Int}, n::Int)
    edges = Tuple{Int,Int}[]
    
    # Add edge from start_check to data1
    push!(edges, ordered_edge(n + start_check, data1))
    
    # Add edges along the path
    for idx in 1:(length(path)-1)
        v1, v2 = path[idx], path[idx+1]
        # Determine if nodes are data or check bits based on the alternating pattern
        if idx % 2 == 1  # v1 is data, v2 is check
            push!(edges, ordered_edge(v1, n + v2))
        else  # v1 is check, v2 is data
            push!(edges, ordered_edge(n + v1, v2))
        end
    end
    
    # Add edge from data2 back to start_check
    push!(edges, ordered_edge(data2, n + start_check))
    
    return edges
end

"""
    get_loops_for_check_bit(H, check_bit_index; max_length=20)

Main API function to get loops in the Tanner graph supported on a given check bit.

Arguments:
- `H`: Parity check matrix (m × n) 
- `check_bit_index`: Check bit index (1 to m)
- `max_length`: Maximum loop length to search (default: 20)

Returns:
- Vector of vectors, where each inner vector represents a loop as a list of edges
  Each edge is a tuple (v1, v2) where:
  - Data bits are numbered 1:n  
  - Check bits are numbered (n+1):(n+m)
  - Loop is represented as [(v1,v2), (v2,v3), ..., (vn,v1)]

Example:
```julia
H = [1 1 0; 1 0 1; 0 1 1]  # 3x3 parity check matrix
loops = get_loops_for_check_bit(H, 1)
# Returns: [[(1, 4), (4, 2), (2, 5), (5, 3), (3, 6), (6, 1)]]
```
"""
function get_loops_for_check_bit(H::AbstractMatrix, check_bit_index::Int; max_length::Int=20)
    # Find all TannerLoop objects
    tanner_loops = find_loops_for_check_bit(H, check_bit_index, max_length=max_length)
    
    # Convert to simple edge lists
    return [loop.edges for loop in tanner_loops]
end
# Make the main function available when this file is included
export TannerLoop, find_tanner_loops, find_loops_for_check_bit, get_loops_for_data_bit, get_loops_for_check_bit, print_loop, test_tanner_loops, test_tanner_loops_comprehensive, validate_loop, create_simple_example, demo_tanner_loops
