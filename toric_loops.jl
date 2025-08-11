"""
Toric Code Loop Finding Module

This module provides functions to find loops in toric code factor graphs.
The main function is `find_toric_code_loops` which finds all loops containing
a specified data qubit up to a given maximum order.

Author: Generated for toric code tensor network belief propagation
"""

module ToricLoops

using SparseArrays

# Export the main structures and functions
export ToricCodeLoop, find_toric_code_loops, find_all_toric_loops_comprehensive

"""
Structure representing a loop in the toric code factor graph.

Fields:
- edges: Vector of (v,w) tuples with v < w representing factor graph edges
- data_qubits: Vector of data qubit indices involved (1 to 2L²)
- check_qubits: Vector of check qubit indices involved (2L²+1 to 3L²)
"""
struct ToricCodeLoop
    edges::Vector{Tuple{Int,Int}}     # Edges (v1,v2) with v1 < v2
    data_qubits::Vector{Int}          # Data qubit indices (1 to 2L²)  
    check_qubits::Vector{Int}         # Check qubit indices (2L²+1 to 3L²)
end

"""
    find_toric_code_loops(pcmat, data_qubit, max_loop_order)

Find all loops up to max_loop_order that contain the specified data_qubit in a toric code.

# Arguments
- `pcmat::SparseMatrixCSC`: Parity check matrix from toric_code_X_parity_matrix(L)
- `data_qubit::Int`: Data qubit index (1 to 2L²)
- `max_loop_order::Int`: Maximum number of edges in loops to find

# Returns
- `Vector{ToricCodeLoop}`: All loops containing the data_qubit

The returned loops have:
- edges: List of (v,w) pairs with v < w representing factor graph edges
- data_qubits: Indices of data qubits involved (1 to 2L²)  
- check_qubits: Indices of check qubits involved (2L²+1 to 3L²)

# Example
```julia
using ToricLoops
L = 3
pcmat = toric_code_X_parity_matrix(L)
loops = find_toric_code_loops(pcmat, 1, 6)
```
"""
function find_toric_code_loops(pcmat::SparseMatrixCSC, data_qubit::Int, max_loop_order::Int)
    m, n = size(pcmat)  # m = L^2 check qubits, n = 2L^2 data qubits
    
    # Build bipartite factor graph adjacency
    data_to_checks = [Int[] for _ in 1:n]
    check_to_data = [Int[] for _ in 1:m] 
    
    for i in 1:m
        for j in 1:n
            if pcmat[i, j] == 1
                push!(data_to_checks[j], i + n)  # Offset check indices by n
                push!(check_to_data[i], j)
            end
        end
    end
    
    found_loops = ToricCodeLoop[]
    
    # Breadth-first search for cycles containing data_qubit
    function find_cycles_containing_qubit(start_qubit::Int, max_depth::Int)
        # Queue entries: (current_node, path, is_data_node, depth)
        queue = [(start_qubit, [start_qubit], true, 0)]
        visited_paths = Set{Vector{Int}}()
        
        while !isempty(queue)
            current, path, is_data, depth = popfirst!(queue)
            
            if depth >= max_depth
                continue
            end
            
            # Skip if we've seen this path before (avoid duplicates)
            path_key = sort(copy(path))
            if path_key in visited_paths
                continue
            end
            
            if is_data
                # Current node is a data qubit - move to connected checks
                for check_node in data_to_checks[current]
                    if length(path) >= 2 && check_node == path[end-1]
                        continue  # No immediate backtrack
                    end
                    
                    new_path = vcat(path, [check_node])
                    push!(queue, (check_node, new_path, false, depth + 1))
                end
            else
                # Current node is a check qubit - move to connected data qubits
                check_idx = current - n
                for data_node in check_to_data[check_idx]
                    
                    # Check if we've found a cycle back to start
                    if data_node == start_qubit && length(path) >= 4 && length(path) % 2 == 0
                        # Valid cycle found! Build the loop structure
                        cycle_path = vcat(path, [data_node])
                        
                        cycle_edges = Tuple{Int,Int}[]
                        cycle_data = Int[]
                        cycle_checks = Int[]
                        
                        # Separate data and check qubits
                        for node in cycle_path[1:end-1]  # Exclude duplicate end node
                            if node <= n
                                push!(cycle_data, node)
                            else
                                push!(cycle_checks, node)
                            end
                        end
                        
                        # Create edges between consecutive nodes
                        for i in 1:(length(cycle_path)-1)
                            v, w = min(cycle_path[i], cycle_path[i+1]), max(cycle_path[i], cycle_path[i+1])
                            push!(cycle_edges, (v, w))
                        end
                        
                        # Remove duplicates and sort
                        unique_edges = unique(cycle_edges)
                        sort!(unique_edges)
                        unique_data = unique(cycle_data)
                        sort!(unique_data)
                        unique_checks = unique(cycle_checks)
                        sort!(unique_checks)
                        
                        # Only add if within size limit and not duplicate
                        if length(unique_edges) <= max_loop_order
                            loop = ToricCodeLoop(unique_edges, unique_data, unique_checks)
                            
                            # Check for existing loops with same edge set
                            is_duplicate = any(existing_loop -> Set(existing_loop.edges) == Set(unique_edges), found_loops)
                            
                            if !is_duplicate
                                push!(found_loops, loop)
                                push!(visited_paths, path_key)
                            end
                        end
                        
                    elseif data_node != start_qubit && length(path) < max_depth
                        # Continue exploring
                        if length(path) >= 2 && data_node == path[end-1]
                            continue  # No immediate backtrack
                        end
                        
                        new_path = vcat(path, [data_node])
                        push!(queue, (data_node, new_path, true, depth + 1))
                    end
                end
            end
        end
    end
    
    # Find all cycles containing our data qubit
    find_cycles_containing_qubit(data_qubit, max_loop_order)
    
    return found_loops
end

"""
    find_all_toric_loops_comprehensive(L, max_loop_order)

Find all loops for every data qubit in the toric code.

# Arguments
- `L::Int`: Linear size of toric code
- `max_loop_order::Int`: Maximum number of edges per loop

# Returns
- `Vector{Vector{ToricCodeLoop}}`: loops[i] contains all loops for data qubit i

# Note
This function requires `toric_code_X_parity_matrix(L)` to be defined in the calling scope.
"""
function find_all_toric_loops_comprehensive(L::Int, max_loop_order::Int)
    # This assumes toric_code_X_parity_matrix is available in the calling scope
    pcmat = toric_code_X_parity_matrix(L)
    n = 2 * L^2
    
    all_loops = Vector{ToricCodeLoop}[]
    
    # Note: Removed @showprogress to avoid dependency - you can add it back if needed
    for data_qubit in 1:n
        loops = find_toric_code_loops(pcmat, data_qubit, max_loop_order)
        push!(all_loops, loops)
    end
    
    return all_loops
end

"""
    get_loop_for_existing_format(loop::ToricCodeLoop, n::Int)

Convert ToricCodeLoop to the format expected by existing tensor network functions.

# Arguments
- `loop::ToricCodeLoop`: The loop to convert
- `n::Int`: Number of data qubits (2L²)

# Returns
- `edges`: Vector of tuples for the loop edges
- `data_bits`: Vector of data qubit indices  
- `check_bits`: Vector of check qubit indices (without offset)
"""
function get_loop_for_existing_format(loop::ToricCodeLoop, n::Int)
    edges = loop.edges
    data_bits = loop.data_qubits
    check_bits = [c - n for c in loop.check_qubits]  # Remove offset for compatibility
    
    return edges, data_bits, check_bits
end

# Helper function for compatibility with existing code structure
"""
    convert_to_tanner_loop_format(loops::Vector{ToricCodeLoop}, n::Int)

Convert vector of ToricCodeLoops to format compatible with existing tannerloopslist.

# Arguments
- `loops::Vector{ToricCodeLoop}`: Loops to convert
- `n::Int`: Number of data qubits (2L²)

# Returns
- Vector of named tuples with fields: edges, data_bits, check_bits
"""
function convert_to_tanner_loop_format(loops::Vector{ToricCodeLoop}, n::Int)
    converted_loops = []
    
    for loop in loops
        edges, data_bits, check_bits = get_loop_for_existing_format(loop, n)
        
        # Create named tuple compatible with existing code
        tanner_loop = (
            edges = edges,
            data_bits = data_bits, 
            check_bits = check_bits
        )
        
        push!(converted_loops, tanner_loop)
    end
    
    return converted_loops
end

end  # End of module ToricLoops
