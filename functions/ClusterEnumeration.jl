"""
ClusterEnumeration.jl

Julia implementation for enumerating connected clusters of loops.
"""

include("LoopEnumeration.jl")

struct Cluster
    loop_ids::Vector{Int}
    multiplicities::Dict{Int, Int}
    weight::Int
    total_loops::Int
end

struct ClusterEnumerator
    adj_matrix::Matrix{Int}
    loop_enumerator::LoopEnumerator
    
    function ClusterEnumerator(adj_matrix::Matrix{Int})
        loop_enum = LoopEnumerator(adj_matrix)
        new(adj_matrix, loop_enum)
    end
end

function enumerate_connected_clusters(enumerator::ClusterEnumerator, site::Int, 
                                    max_weight::Int, max_loop_weight::Int = max_weight)
    """
    Enumerate all connected clusters of total weight <= max_weight supported on site.
    """
    
    println("Step 1: Enumerating all loops up to weight $max_loop_weight...")
    
    # Find all possible loops in the graph
    all_loops = find_all_loops_in_graph(enumerator, max_loop_weight)
    println("Found $(length(all_loops)) total loops")
    
    # Filter loops supported on target site
    supported_loops = [loop for loop in all_loops if site in loop.vertices]
    println("Found $(length(supported_loops)) loops supported on site $site")
    
    # Build interaction graph
    println("Step 2: Building interaction graph...")
    interaction_graph = build_interaction_graph(all_loops)
    
    # Enumerate clusters
    println("Step 3: Enumerating clusters...")
    all_clusters = enumerate_clusters_up_to_weight(all_loops, supported_loops, 
                                                  site, max_weight, interaction_graph)
    println("Found $(length(all_clusters)) total clusters before connectivity filtering")
    
    # Filter for connected clusters
    println("Step 4: Filtering for connected clusters...")
    connected_clusters = filter_connected_clusters(all_clusters, interaction_graph)
    println("Found $(length(connected_clusters)) connected clusters")
    
    # Remove redundancies
    println("Step 5: Removing redundancies...")
    unique_clusters = remove_cluster_redundancies(connected_clusters)
    println("Found $(length(unique_clusters)) unique connected clusters")
    
    return unique_clusters
end

function find_all_loops_in_graph(enumerator::ClusterEnumerator, max_weight::Int)
    all_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Vector{Int}}, Int}}()
    loop_id = 0
    
    for vertex in 1:enumerator.loop_enumerator.n_vertices
        vertex_loops = find_loops_supported_on_vertex(enumerator.loop_enumerator, vertex, max_weight)
        
        for loop in vertex_loops
            canonical = canonical_loop_representation(loop)
            
            if !(canonical in seen_loops)
                push!(seen_loops, canonical)
                # Create new loop with ID
                new_loop = Loop(sort(loop.vertices), sort(loop.edges), loop.weight)
                push!(all_loops, new_loop)
                loop_id += 1
            end
        end
    end
    
    return all_loops
end

function canonical_loop_representation(loop::Loop)
    vertices = sort(loop.vertices)
    edges = sort([sort([u, v]) for (u, v) in loop.edges])
    return (vertices, edges, loop.weight)
end

function build_interaction_graph(loops::Vector{Loop})
    interaction_graph = Dict{Int, Vector{Int}}()
    
    for i in 1:length(loops)
        interaction_graph[i] = Int[]
        for j in 1:length(loops)
            if i != j
                vertices1 = Set(loops[i].vertices)
                vertices2 = Set(loops[j].vertices)
                
                if !isempty(intersect(vertices1, vertices2))
                    push!(interaction_graph[i], j)
                end
            end
        end
    end
    
    return interaction_graph
end

function enumerate_clusters_up_to_weight(all_loops::Vector{Loop}, supported_loops::Vector{Loop},
                                       site::Int, max_weight::Int, 
                                       interaction_graph::Dict{Int, Vector{Int}})
    clusters = Cluster[]
    
    # Create mapping from loop to its ID by comparing canonical representations
    supported_loop_ids = Int[]
    for supported_loop in supported_loops
        canonical_supported = canonical_loop_representation(supported_loop)
        for (i, loop) in enumerate(all_loops)
            canonical_loop = canonical_loop_representation(loop)
            if canonical_supported == canonical_loop
                push!(supported_loop_ids, i)
                break
            end
        end
    end
    
    for total_weight in 1:max_weight
        append!(clusters, generate_clusters_of_weight(all_loops, supported_loop_ids, total_weight))
    end
    
    return clusters
end

function generate_clusters_of_weight(all_loops::Vector{Loop}, supported_loop_ids::Vector{Int}, 
                                   target_weight::Int)
    clusters = Cluster[]
    
    function generate_partitions(remaining_weight::Int, start_loop_idx::Int, current_cluster::Vector{Int})
        if remaining_weight == 0
            if any(loop_id in supported_loop_ids for loop_id in current_cluster)
                multiplicities = Dict{Int, Int}()
                for loop_id in current_cluster
                    multiplicities[loop_id] = get(multiplicities, loop_id, 0) + 1
                end
                
                cluster = Cluster(
                    collect(keys(multiplicities)),
                    multiplicities,
                    target_weight,
                    length(current_cluster)
                )
                push!(clusters, cluster)
            end
            return
        end
        
        if remaining_weight < 0
            return
        end
        
        for loop_idx in start_loop_idx:length(all_loops)
            loop_weight = all_loops[loop_idx].weight
            
            if loop_weight <= remaining_weight
                push!(current_cluster, loop_idx)
                generate_partitions(remaining_weight - loop_weight, loop_idx, current_cluster)
                pop!(current_cluster)
            end
        end
    end
    
    generate_partitions(target_weight, 1, Int[])
    return clusters
end

function filter_connected_clusters(clusters::Vector{Cluster}, 
                                 interaction_graph::Dict{Int, Vector{Int}})
    connected_clusters = Cluster[]
    
    for cluster in clusters
        if is_cluster_connected(cluster, interaction_graph)
            push!(connected_clusters, cluster)
        end
    end
    
    return connected_clusters
end

function is_cluster_connected(cluster::Cluster, interaction_graph::Dict{Int, Vector{Int}})
    loop_ids = cluster.loop_ids
    multiplicities = cluster.multiplicities
    
    if length(loop_ids) == 1
        return true
    end
    
    # Create vertices for each loop instance
    vertices = Tuple{Int, Int}[]
    vertex_id = 1
    
    for loop_id in loop_ids
        mult = multiplicities[loop_id]
        for _ in 1:mult
            push!(vertices, (loop_id, vertex_id))
            vertex_id += 1
        end
    end
    
    if length(vertices) <= 1
        return true
    end
    
    # Build adjacency list for cluster interaction graph
    cluster_adj = Dict{Int, Vector{Int}}()
    
    for i in 1:length(vertices)
        cluster_adj[vertices[i][2]] = Int[]
        for j in 1:length(vertices)
            if i != j
                loop_id1, v1 = vertices[i]
                loop_id2, v2 = vertices[j]
                
                if loop_id2 in get(interaction_graph, loop_id1, Int[])
                    push!(cluster_adj[v1], v2)
                end
            end
        end
    end
    
    # Check connectivity using BFS
    visited = Set{Int}()
    start_vertex = vertices[1][2]
    queue = [start_vertex]
    push!(visited, start_vertex)
    
    while !isempty(queue)
        current = popfirst!(queue)
        for neighbor in cluster_adj[current]
            if !(neighbor in visited)
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    
    return length(visited) == length(vertices)
end

function remove_cluster_redundancies(clusters::Vector{Cluster})
    unique_clusters = Cluster[]
    seen_signatures = Set{Tuple}()
    
    for cluster in clusters
        signature = canonical_cluster_signature(cluster)
        
        if !(signature in seen_signatures)
            push!(seen_signatures, signature)
            push!(unique_clusters, cluster)
        end
    end
    
    return unique_clusters
end

function canonical_cluster_signature(cluster::Cluster)
    items = [(loop_id, cluster.multiplicities[loop_id]) for loop_id in sort(cluster.loop_ids)]
    return (tuple(items...), cluster.weight)
end

function ursell_function(cluster::Cluster, all_loops::Vector{Loop})
    """
    Compute the Ursell function φ(W) for a given cluster W.
    
    According to the polymer proof:
    - φ(W) = 0 if the cluster is disconnected
    - φ(W) = (1/W!) * Σ_{C connected subgraphs of G_W} (-1)^{|E(C)|} if connected
    
    Args:
        cluster: The cluster W = {(l_i, η_i)} 
        all_loops: Vector of all loops to build interaction graph
    
    Returns:
        Float64: The Ursell function value
    """
    
    # Check if cluster is connected using existing function
    interaction_graph = build_interaction_graph(all_loops)
    
    if !is_cluster_connected(cluster, interaction_graph)
        return 0.0  # Disconnected clusters have Ursell function = 0
    end
    
    # For connected clusters, compute the sum over connected subgraphs
    cluster_interaction_graph = build_cluster_interaction_graph(cluster, all_loops)
    
    # Compute 1/W! where W! = ∏_i η_i!
    w_factorial = 1.0
    for multiplicity in values(cluster.multiplicities)
        w_factorial *= factorial(multiplicity)
    end
    
    # Sum over all connected subgraphs of the cluster interaction graph
    subgraph_sum = sum_over_connected_subgraphs(cluster_interaction_graph)
    
    return subgraph_sum / w_factorial
end

function build_cluster_interaction_graph(cluster::Cluster, all_loops::Vector{Loop})
    """
    Build the interaction graph G_W for a cluster W.
    
    Vertices are loop instances (accounting for multiplicities).
    Edges connect loops that are incompatible (share vertices) or identical.
    """
    
    # Create vertex list: each loop appears η_i times
    vertices = Tuple{Int, Int}[]  # (loop_id, instance_id)
    vertex_id = 1
    
    for loop_id in cluster.loop_ids
        multiplicity = cluster.multiplicities[loop_id]
        for instance in 1:multiplicity
            push!(vertices, (loop_id, vertex_id))
            vertex_id += 1
        end
    end
    
    # Build adjacency list for cluster interaction graph
    adj_list = Dict{Int, Vector{Int}}()
    
    for i in 1:length(vertices)
        adj_list[i] = Int[]
        for j in 1:length(vertices)
            if i != j
                loop_id1, _ = vertices[i]
                loop_id2, _ = vertices[j]
                
                # Connect if loops are identical or incompatible (share vertices)
                if loop_id1 == loop_id2 || loops_are_incompatible(all_loops[loop_id1], all_loops[loop_id2])
                    push!(adj_list[i], j)
                end
            end
        end
    end
    
    return adj_list
end

function loops_are_incompatible(loop1::Loop, loop2::Loop)
    """Check if two loops are incompatible (share at least one vertex)."""
    vertices1 = Set(loop1.vertices)
    vertices2 = Set(loop2.vertices)
    return !isempty(intersect(vertices1, vertices2))
end

function sum_over_connected_subgraphs(adj_list::Dict{Int, Vector{Int}})
    """
    Compute Σ_{C connected subgraphs} (-1)^{|E(C)|}.
    
    This uses the inclusion-exclusion principle over all connected subgraphs.
    """
    n_vertices = length(adj_list)
    
    if n_vertices == 0
        return 1.0  # Empty graph
    end
    
    if n_vertices == 1
        return 1.0  # Single vertex, no edges
    end
    
    total_sum = 0.0
    
    # Iterate over all non-empty subsets of vertices
    for subset_bits in 1:(2^n_vertices - 1)
        subset = Int[]
        for i in 1:n_vertices
            if (subset_bits >> (i-1)) & 1 == 1
                push!(subset, i)
            end
        end
        
        # Check if this subset forms a connected subgraph
        if is_subgraph_connected(subset, adj_list)
            # Count edges in this subgraph
            edge_count = count_edges_in_subgraph(subset, adj_list)
            
            # Add (-1)^{|E(C)|} to the sum
            total_sum += (-1.0)^edge_count
        end
    end
    
    return total_sum
end

function is_subgraph_connected(vertices::Vector{Int}, adj_list::Dict{Int, Vector{Int}})
    """Check if a subset of vertices forms a connected subgraph."""
    if length(vertices) <= 1
        return true
    end
    
    vertex_set = Set(vertices)
    visited = Set{Int}()
    queue = [vertices[1]]
    push!(visited, vertices[1])
    
    while !isempty(queue)
        current = popfirst!(queue)
        for neighbor in adj_list[current]
            if neighbor in vertex_set && !(neighbor in visited)
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    
    return length(visited) == length(vertices)
end

function count_edges_in_subgraph(vertices::Vector{Int}, adj_list::Dict{Int, Vector{Int}})
    """Count the number of edges in a subgraph."""
    vertex_set = Set(vertices)
    edge_count = 0
    
    for v in vertices
        for neighbor in adj_list[v]
            if neighbor in vertex_set && v < neighbor  # Avoid double counting
                edge_count += 1
            end
        end
    end
    
    return edge_count
end

# Test function for Ursell function
function test_ursell_function()
    """Test that disconnected clusters have Ursell function = 0."""
    println("Testing Ursell function...")
    
    # Create a simple 4x4 lattice for testing
    L = 4
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = ClusterEnumerator(adj_matrix)
    
    # Find some loops
    all_loops = find_all_loops_in_graph(enumerator, 6)
    println("Found $(length(all_loops)) loops for testing")
    
    if length(all_loops) >= 2
        # Create a disconnected cluster (two loops that don't share vertices)
        loop1 = all_loops[1]
        loop2 = nothing
        
        # Find a loop that doesn't share vertices with loop1
        for loop in all_loops[2:end]
            if !loops_are_incompatible(loop1, loop)
                loop2 = loop
                break
            end
        end
        
        if loop2 !== nothing
            # Create disconnected cluster
            disconnected_cluster = Cluster(
                [1, 2],  # loop IDs
                Dict(1 => 1, 2 => 1),  # multiplicities
                loop1.weight + loop2.weight,  # total weight
                2  # total loops
            )
            
            ursell_val = ursell_function(disconnected_cluster, [loop1, loop2])
            println("Disconnected cluster Ursell function: $ursell_val")
            
            if abs(ursell_val) < 1e-10
                println("✅ Test PASSED: Disconnected cluster has Ursell function ≈ 0")
            else
                println("❌ Test FAILED: Expected Ursell function ≈ 0, got $ursell_val")
            end
        else
            println("⚠️  Could not find two non-overlapping loops for disconnected test")
        end
        
        # Test a simple connected cluster (single loop)
        single_cluster = Cluster(
            [1],  # single loop
            Dict(1 => 1),  # multiplicity 1
            loop1.weight,  # weight
            1  # one loop
        )
        
        ursell_val = ursell_function(single_cluster, [loop1])
        println("Single loop cluster Ursell function: $ursell_val")
        println("Expected: 1.0 (single vertex graph has sum = 1)")
        
        if abs(ursell_val - 1.0) < 1e-10
            println("✅ Test PASSED: Single loop has Ursell function = 1")
        else
            println("❌ Test FAILED: Expected Ursell function = 1, got $ursell_val")
        end
    else
        println("⚠️  Not enough loops found for testing")
    end
end

# No module exports needed for include-style usage