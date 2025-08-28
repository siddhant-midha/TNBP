"""
ClusterEnumeration.jl

Julia implementation for enumerating connected clusters of loops.
    Replaced with simpler Ursell function for safety
"""

include("LoopEnumeration.jl")
using ProgressMeter
using Serialization
using Dates

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
    Uses DFS starting from loops supported on the target site.
    """
    
    println("Step 1: Enumerating all loops up to weight $max_loop_weight...")
    
    # Find all possible loops in the graph
    # all_loops = find_all_loops_in_graph(enumerator, max_loop_weight)
    all_loops = find_all_loops_in_graph(enumerator, max_loop_weight)
    println("Found $(length(all_loops)) total loops")
    
    # Build interaction graph once
    println("Step 2: Building interaction graph...")
    interaction_graph = build_interaction_graph_optimized(all_loops)
    
    # Find loops supported on target site (these are our starting points)
    supported_loop_ids = Int[]
    for (i, loop) in enumerate(all_loops)
        if site in loop.vertices
            push!(supported_loop_ids, i)
        end
    end
    println("Found $(length(supported_loop_ids)) loops supported on site $site")
    
    if isempty(supported_loop_ids)
        println("No loops found on site $site - returning empty cluster list")
        return Cluster[]
    end
    
    # DFS cluster enumeration starting from supported loops
    println("Step 3: DFS cluster enumeration from supported loops...")
    all_clusters = dfs_enumerate_clusters_from_supported(all_loops, supported_loop_ids, 
                                                        max_weight, interaction_graph)
    println("Found $(length(all_clusters)) connected clusters")
    
    # Remove redundancies
    println("Step 4: Removing redundancies...")
    unique_clusters = remove_cluster_redundancies(all_clusters)
    println("Found $(length(unique_clusters)) unique connected clusters")
    
    return unique_clusters
end

function find_all_loops_in_graph(enumerator::ClusterEnumerator, max_weight::Int)
    """Find all loops by enumerating from each vertex but avoiding redundant work."""
    all_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    n_vertices = enumerator.loop_enumerator.n_vertices
    println("  Enumerating loops from $n_vertices vertices...")
    
    # Progress bar for loop enumeration
    println("  Starting loop enumeration...")
    flush(stdout)
    progress = Progress(n_vertices, dt=0.1, desc="Finding loops: ", color=:blue, barlen=50)
    
    # Only enumerate from the smallest vertex in each potential loop to avoid redundancy
    for vertex in 1:n_vertices
        vertex_loops = find_loops_supported_on_vertex_optimized(enumerator.loop_enumerator, vertex, max_weight, seen_loops)
        # vertex_loops = find_loops_supported_on_vertex(enumerator.loop_enumerator, vertex, max_weight)
        
        for loop in vertex_loops
            # Only keep loops where 'vertex' is the smallest vertex (canonical representative)
            if vertex == minimum(loop.vertices)
                canonical = canonical_loop_representation(loop)
                
                if !(canonical in seen_loops)
                    push!(seen_loops, canonical)
                    new_loop = Loop(sort(loop.vertices), sort(loop.edges), loop.weight)
                    push!(all_loops, new_loop)
                end
            end
        end
        
        next!(progress, showvalues = [("Vertex", "$vertex/$n_vertices"), ("Loops found", "$(length(all_loops))")])
    end
    
    return all_loops
end

# Import canonical_loop_representation from LoopEnumeration.jl (no need to redefine)

function build_interaction_graph(loops::Vector{Loop})
    interaction_graph = Dict{Int, Vector{Int}}()
    n_loops = length(loops)
    
    println("  Building interaction graph for $n_loops loops...")
    progress = Progress(n_loops, dt=1.0, desc="Building graph: ")
    
    for i in 1:n_loops
        interaction_graph[i] = Int[]
        for j in 1:n_loops
            if i != j
                vertices1 = Set(loops[i].vertices)
                vertices2 = Set(loops[j].vertices)
                
                if !isempty(intersect(vertices1, vertices2))
                    push!(interaction_graph[i], j)
                end
            end
        end
        next!(progress)
    end
    
    return interaction_graph
end

function build_interaction_graph_optimized(loops::Vector{Loop})
    """Build interaction graph with optimizations for speed."""
    interaction_graph = Dict{Int, Vector{Int}}()
    n_loops = length(loops)
    
    println("  Building optimized interaction graph for $n_loops loops...")
    flush(stdout)
    progress = Progress(n_loops, dt=0.1, desc="Building graph: ", color=:green, barlen=50)
    
    # Optimization 1: Pre-compute vertex sets once
    vertex_sets = [Set(loop.vertices) for loop in loops]
    
    # Optimization 2: Build vertex-to-loops mapping for faster lookup
    vertex_to_loops = Dict{Int, Vector{Int}}()
    for (i, loop) in enumerate(loops)
        for vertex in loop.vertices
            if !haskey(vertex_to_loops, vertex)
                vertex_to_loops[vertex] = Int[]
            end
            push!(vertex_to_loops[vertex], i)
        end
    end
    
    for i in 1:n_loops
        interaction_graph[i] = Int[]
        
        # Optimization 3: Only check loops that share at least one vertex
        candidate_loops = Set{Int}()
        for vertex in loops[i].vertices
            for j in get(vertex_to_loops, vertex, Int[])
                if i != j
                    push!(candidate_loops, j)
                end
            end
        end
        
        # Now check intersection only with candidates
        for j in candidate_loops
            if !isempty(intersect(vertex_sets[i], vertex_sets[j]))
                push!(interaction_graph[i], j)
            end
        end
        
        next!(progress)
    end
    
    return interaction_graph
end

function dfs_enumerate_clusters_from_supported(all_loops::Vector{Loop}, supported_loop_ids::Vector{Int}, 
                                               max_weight::Int, interaction_graph::Dict{Int, Vector{Int}})
    """
    Enumerate connected clusters using DFS starting from loops supported on target site.
    Connectivity is guaranteed by growing through the interaction graph.
    """
    clusters = Cluster[]
    seen_clusters = Set{Tuple}()
    cluster_count = 0
    
    println("  Starting DFS cluster enumeration...")
    println("  Supported loops: $(length(supported_loop_ids)), Max weight: $max_weight")
    
    # Progress tracking
    last_report_time = time()
    last_cluster_count = 0
    
    # DFS to grow clusters starting from each supported loop
    function dfs_grow_cluster(current_cluster::Vector{Int}, current_weight::Int, 
                             has_supported::Bool)
        
        # If we've found a valid cluster (has supported loop), record it
        if has_supported && current_weight >= 1
            # Create cluster with multiplicities
            multiplicities = Dict{Int, Int}()
            for loop_id in current_cluster
                multiplicities[loop_id] = get(multiplicities, loop_id, 0) + 1
            end
            
            cluster = Cluster(
                collect(keys(multiplicities)),
                multiplicities,
                current_weight,
                length(current_cluster)
            )
            
            # Avoid duplicates using canonical signature
            signature = canonical_cluster_signature(cluster)
            if !(signature in seen_clusters)
                push!(seen_clusters, signature)
                push!(clusters, cluster)
                cluster_count += 1
                
                # Progress reporting every 2 seconds
                current_time = time()
                if current_time - last_report_time >= 2.0
                    new_clusters = cluster_count - last_cluster_count
                    println("    Found $cluster_count clusters (+$new_clusters in last 2s)")
                    last_report_time = current_time
                    last_cluster_count = cluster_count
                end
            end
        end
        
        # Stop if we've reached max weight
        if current_weight >= max_weight
            return
        end
        
        # Find candidate loops to add (adjacent loops or multiplicities)
        candidate_loops = Set{Int}()
        
        if isempty(current_cluster)
            # Start with supported loops only
            for loop_id in supported_loop_ids
                if all_loops[loop_id].weight <= max_weight - current_weight
                    push!(candidate_loops, loop_id)
                end
            end
        else
            # Add loops connected to current cluster via interaction graph
            for loop_id in current_cluster
                # Add connected loops (touching loops)
                for neighbor_id in get(interaction_graph, loop_id, Int[])
                    if all_loops[neighbor_id].weight <= max_weight - current_weight
                        push!(candidate_loops, neighbor_id)
                    end
                end
                # Allow multiplicity increases (same loop added again)
                if all_loops[loop_id].weight <= max_weight - current_weight
                    push!(candidate_loops, loop_id)
                end
            end
        end
        
        # Try each candidate loop
        for loop_id in candidate_loops
            loop_weight = all_loops[loop_id].weight
            new_weight = current_weight + loop_weight
            
            if new_weight <= max_weight
                new_cluster = copy(current_cluster)
                push!(new_cluster, loop_id)
                new_has_supported = has_supported || (loop_id in supported_loop_ids)
                
                # Continue DFS (connectivity guaranteed by interaction graph)
                dfs_grow_cluster(new_cluster, new_weight, new_has_supported)
            end
        end
    end
    
    # Start DFS with empty cluster
    dfs_grow_cluster(Int[], 0, false)
    
    println("  DFS enumeration completed: $cluster_count total clusters found")
    return clusters
end

function generate_clusters_of_weight(all_loops::Vector{Loop}, supported_loop_ids::Vector{Int}, 
                                   target_weight::Int)
    """Generate clusters efficiently by building connected components incrementally."""
    clusters = Cluster[]
    
    # Build interaction graph once for efficiency
    interaction_graph = build_interaction_graph(all_loops)
    
    function generate_connected_partitions(remaining_weight::Int, current_cluster::Vector{Int}, 
                                         used_loop_ids::Set{Int}, must_include_supported::Bool)
        if remaining_weight == 0
            if must_include_supported && any(loop_id in supported_loop_ids for loop_id in current_cluster)
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
                
                # Only add if connected (check connectivity during construction)
                if is_cluster_connected(cluster, interaction_graph)
                    push!(clusters, cluster)
                end
            end
            return
        end
        
        if remaining_weight < 0
            return
        end
        
        # Find extendable loops (connected to current cluster or starting new)
        candidate_loops = Int[]
        
        if isempty(current_cluster)
            # Start with any loop
            for loop_idx in 1:length(all_loops)
                if all_loops[loop_idx].weight <= remaining_weight
                    push!(candidate_loops, loop_idx)
                end
            end
        else
            # Only add loops connected to current cluster or same loops (multiplicities)
            connected_loops = Set{Int}()
            for loop_id in current_cluster
                # Add connected loops
                for connected_id in get(interaction_graph, loop_id, Int[])
                    if !(connected_id in used_loop_ids) && all_loops[connected_id].weight <= remaining_weight
                        push!(connected_loops, connected_id)
                    end
                end
                # Allow multiplicity increases
                if all_loops[loop_id].weight <= remaining_weight
                    push!(connected_loops, loop_id)
                end
            end
            candidate_loops = collect(connected_loops)
        end
        
        # Try each candidate
        for loop_idx in candidate_loops
            loop_weight = all_loops[loop_idx].weight
            if loop_weight <= remaining_weight
                new_cluster = copy(current_cluster)
                push!(new_cluster, loop_idx)
                new_used = copy(used_loop_ids)
                if !(loop_idx in current_cluster)  # Only mark as used if it's truly new
                    push!(new_used, loop_idx)
                end
                new_must_include = must_include_supported || (loop_idx in supported_loop_ids)
                
                generate_connected_partitions(remaining_weight - loop_weight, new_cluster, 
                                            new_used, new_must_include)
            end
        end
    end
    
    generate_connected_partitions(target_weight, Int[], Set{Int}(), false)
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

function create_coordinate_canonical_form_for_cluster(loop::Loop, L::Int = 11)
    """
    Create a true canonical form based on coordinate geometry for cluster identification.
    This handles periodic boundary conditions properly for translation symmetry tests.
    """
    # Convert vertices to coordinates
    function site_to_coords(site::Int)
        i = div(site - 1, L) + 1
        j = mod(site - 1, L) + 1
        return (i, j)
    end
    
    coords = [site_to_coords(v) for v in sort(loop.vertices)]
    
    # Check for periodic wrapping by examining coordinate spans
    i_values = [coord[1] for coord in coords]
    j_values = [coord[2] for coord in coords]
    
    i_span = maximum(i_values) - minimum(i_values)
    j_span = maximum(j_values) - minimum(j_values)
    
    # If span > L/2, the pattern likely wraps around periodic boundaries
    # We need to "unwrap" it to get the true geometric pattern
    if i_span > L/2 || j_span > L/2
        # Unwrap coordinates by shifting wrapped vertices
        unwrapped_coords = []
        
        for coord in coords
            i, j = coord
            # Unwrap i-coordinate if needed
            if i_span > L/2
                # Find vertices that are "far" from the minimum i
                min_i = minimum(i_values)
                if i - min_i > L/2
                    i = i - L  # Shift back to represent wrapping
                end
            end
            
            # Unwrap j-coordinate if needed  
            if j_span > L/2
                # Find vertices that are "far" from the minimum j
                min_j = minimum(j_values)
                if j - min_j > L/2
                    j = j - L  # Shift back to represent wrapping
                end
            end
            
            push!(unwrapped_coords, (i, j))
        end
        
        coords = unwrapped_coords
    end
    
    # Normalize coordinates to start from (0, 0)
    min_i = minimum(coord[1] for coord in coords)
    min_j = minimum(coord[2] for coord in coords)
    
    normalized_coords = [(coord[1] - min_i, coord[2] - min_j) for coord in coords]
    
    # Sort for canonical ordering
    normalized_coords = sort(normalized_coords)
    
    return (normalized_coords, loop.weight)
end

function translation_aware_cluster_signature(cluster::Cluster, all_loops::Vector{Loop}, L::Int = 11)
    """
    Create cluster signature using translation-aware canonical loop representations.
    This properly identifies clusters that are translations of each other.
    """
    canonical_loop_sigs = []
    for loop_id in cluster.loop_ids
        loop = all_loops[loop_id]
        # Create proper coordinate-based canonical form
        canonical_loop = create_coordinate_canonical_form_for_cluster(loop, L)
        multiplicity = cluster.multiplicities[loop_id]
        push!(canonical_loop_sigs, (canonical_loop, multiplicity))
    end
    
    # Sort by canonical loop representation for consistent ordering
    sort!(canonical_loop_sigs)
    
    return (tuple(canonical_loop_sigs...), cluster.weight)
end
# COMMENTING OUT FOR SAFETY
# function ursell_function(cluster::Cluster, all_loops::Vector{Loop})
#     """
#     Compute the Ursell function φ(W) for a given cluster W.
    
#     According to the polymer proof:
#     - φ(W) = 0 if the cluster is disconnected
#     - φ(W) = (1/W!) * Σ_{C connected subgraphs of G_W} (-1)^{|E(C)|} if connected
    
#     Args:
#         cluster: The cluster W = {(l_i, η_i)} 
#         all_loops: Vector of all loops to build interaction graph
    
#     Returns:
#         Float64: The Ursell function value
#     """
    
#     # Check if cluster is connected using existing function
#     interaction_graph = build_interaction_graph(all_loops)
    
#     if !is_cluster_connected(cluster, interaction_graph)
#         return 0.0  # Disconnected clusters have Ursell function = 0
#     end
    
#     # For connected clusters, compute the sum over connected subgraphs
#     cluster_interaction_graph = build_cluster_interaction_graph(cluster, all_loops)
    
#     # Compute 1/W! where W! = ∏_i η_i!
#     w_factorial = 1.0
#     for multiplicity in values(cluster.multiplicities)
#         w_factorial *= factorial(multiplicity)
#     end
    
#     # Sum over all connected subgraphs of the cluster interaction graph
#     subgraph_sum = sum_over_connected_subgraphs(cluster_interaction_graph)
    
#     return subgraph_sum / w_factorial
# end

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

# Cluster enumeration and saving/loading functions

struct ClusterEnumerationData
    """Data structure to save complete cluster enumeration results."""
    clusters_by_site::Dict{Int, Vector{Cluster}}
    all_loops::Vector{Loop}
    adj_matrix::Matrix{Int}
    max_weight::Int
    lattice_size::Int
    enumeration_time::Float64
    timestamp::String
    total_sites::Int
end

function enumerate_all_clusters_all_sites(enumerator::ClusterEnumerator, max_weight::Int, 
                                         max_loop_weight::Int = max_weight; 
                                         prefix::String = "")
    """
    Enumerate all connected clusters on all sites and save results.
    
    Args:
        enumerator: ClusterEnumerator with the lattice graph
        max_weight: Maximum cluster weight to enumerate
        max_loop_weight: Maximum individual loop weight (defaults to max_weight)
        prefix: Optional prefix for filename
    
    Returns:
        ClusterEnumerationData: Complete enumeration results
    """
    
    println("🚀 Starting complete cluster enumeration on all sites")
    println("="^70)
    
    # Get lattice info
    n_sites = enumerator.loop_enumerator.n_vertices
    L = enumerator.loop_enumerator.L
    
    println("Parameters:")
    println("  Lattice size: $(L)x$(L) ($n_sites sites)")
    println("  Max cluster weight: $max_weight")
    println("  Max loop weight: $max_loop_weight")
    println("  Will save to: saved_clusters/")
    
    start_time = time()
    
    # Find all loops once (shared across all sites)
    println("\n📊 Step 1: Enumerating all loops...")
    println("   This may take a while for larger lattices...")
    flush(stdout)
    all_loops = find_all_loops_in_graph(enumerator, max_loop_weight)
    println("Found $(length(all_loops)) total loops")
    
    # Build interaction graph once 
    println("\n🔗 Step 2: Building interaction graph...")
    flush(stdout)
    interaction_graph = build_interaction_graph_optimized(all_loops)
    
    # Enumerate clusters for each site
    println("\n🌐 Step 3: Enumerating clusters for each site...")
    flush(stdout)
    clusters_by_site = Dict{Int, Vector{Cluster}}()
    
    println("  Processing $(n_sites) sites...")
    flush(stdout)
    progress = Progress(n_sites, dt=0.1, desc="Sites processed: ", color=:cyan, barlen=50)
    
    for site in 1:n_sites
        # Find loops supported on this site
        supported_loop_ids = Int[]
        for (i, loop) in enumerate(all_loops)
            if site in loop.vertices
                push!(supported_loop_ids, i)
            end
        end
        
        if !isempty(supported_loop_ids)
            # Enumerate clusters for this site
            site_clusters = dfs_enumerate_clusters_from_supported(all_loops, supported_loop_ids, 
                                                                max_weight, interaction_graph)
            # Remove redundancies
            unique_clusters = remove_cluster_redundancies(site_clusters)
            clusters_by_site[site] = unique_clusters
        else
            clusters_by_site[site] = Cluster[]
        end
        
        next!(progress)
    end
    
    end_time = time()
    enumeration_time = end_time - start_time
    
    # Create data structure
    data = ClusterEnumerationData(
        clusters_by_site,
        all_loops,
        enumerator.adj_matrix,
        max_weight,
        L,
        enumeration_time,
        string(now()),
        n_sites
    )
    
    # Save results
    save_cluster_enumeration(data, prefix)
    
    # Print summary
    total_clusters = sum(length(clusters) for clusters in values(clusters_by_site))
    println("\n✅ ENUMERATION COMPLETE!")
    println("="^50)
    println("Total enumeration time: $(round(enumeration_time, digits=2)) seconds")
    println("Total clusters found: $total_clusters")
    println("Average clusters per site: $(round(total_clusters / n_sites, digits=2))")
    
    # Show distribution
    cluster_counts = [length(clusters) for clusters in values(clusters_by_site)]
    println("Cluster count distribution:")
    println("  Min: $(minimum(cluster_counts))")
    println("  Max: $(maximum(cluster_counts))")
    println("  Mean: $(round(sum(cluster_counts) / length(cluster_counts), digits=2))")
    
    return data
end

function save_cluster_enumeration(data::ClusterEnumerationData, prefix::String = "")
    """Save cluster enumeration data to file."""
    
    # Create directory if it doesn't exist
    save_dir = "saved_clusters"
    if !isdir(save_dir)
        mkpath(save_dir)
        println("📁 Created directory: $save_dir")
    end
    
    # Generate filename without timestamp
    base_name = "clusters_L$(data.lattice_size)_w$(data.max_weight)"
    if !isempty(prefix)
        base_name = "$(prefix)_$(base_name)"
    end
    
    filepath = joinpath(save_dir, "$(base_name).jld2")
    
    # Save using Julia's serialization
    println("💾 Saving enumeration data to: $(filepath)")
    
    # Create summary for quick inspection
    summary = Dict(
        "lattice_size" => data.lattice_size,
        "total_sites" => data.total_sites,
        "max_weight" => data.max_weight,
        "total_loops" => length(data.all_loops),
        "total_clusters" => sum(length(clusters) for clusters in values(data.clusters_by_site)),
        "enumeration_time_seconds" => data.enumeration_time,
        "timestamp" => data.timestamp
    )
    
    # Save both data and summary
    save_data = Dict(
        "data" => data,
        "summary" => summary
    )
    
    open(filepath, "w") do io
        serialize(io, save_data)
    end
    
    println("✅ Saved successfully!")
    println("   Summary: $summary")
    
    return filepath
end

function load_cluster_enumeration(filepath::String)
    """
    Load saved cluster enumeration data.
    
    Args:
        filepath: Path to saved .jld2 file
        
    Returns:
        ClusterEnumerationData: Loaded enumeration results
    """
    
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    
    println("📖 Loading cluster enumeration from: $filepath")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    
    data = loaded_data["data"]
    summary = loaded_data["summary"]
    
    println("✅ Loaded successfully!")
    println("   Summary: $summary")
    
    return data
end

function list_saved_cluster_files(save_dir::String = "saved_clusters")
    """List all saved cluster enumeration files with their summaries."""
    
    if !isdir(save_dir)
        println("Directory does not exist: $save_dir")
        return
    end
    
    files = filter(f -> endswith(f, ".jld2"), readdir(save_dir))
    
    if isempty(files)
        println("No saved cluster files found in: $save_dir")
        return
    end
    
    println("📋 Saved cluster enumeration files:")
    println("="^60)
    
    for file in files
        filepath = joinpath(save_dir, file)
        try
            loaded_data = open(filepath, "r") do io
                deserialize(io)
            end
            summary = loaded_data["summary"]
            
            println("📄 $file")
            println("   Lattice: $(summary["lattice_size"])x$(summary["lattice_size"]) ($(summary["total_sites"]) sites)")
            println("   Max weight: $(summary["max_weight"])")
            println("   Clusters: $(summary["total_clusters"]), Loops: $(summary["total_loops"])")
            println("   Time: $(round(summary["enumeration_time_seconds"], digits=2))s")
            println("   Created: $(summary["timestamp"])")
            println()
        catch e
            println("❌ Error reading $file: $e")
        end
    end
end

function get_clusters_for_site(data::ClusterEnumerationData, site::Int)
    """Get all clusters supported on a specific site."""
    return get(data.clusters_by_site, site, Cluster[])
end

function get_clusters_by_weight(data::ClusterEnumerationData, weight::Int)
    """Get all clusters of a specific weight across all sites."""
    result = Cluster[]
    
    for clusters in values(data.clusters_by_site)
        for cluster in clusters
            if cluster.weight == weight
                push!(result, cluster)
            end
        end
    end
    
    return result
end

function analyze_saved_enumeration(data::ClusterEnumerationData)
    """Analyze and display statistics about saved enumeration data."""
    
    println("📊 CLUSTER ENUMERATION ANALYSIS")
    println("="^50)
    
    # Basic info
    println("Lattice: $(data.lattice_size)x$(data.lattice_size) ($(data.total_sites) sites)")
    println("Max weight: $(data.max_weight)")
    println("Total loops: $(length(data.all_loops))")
    println("Enumeration time: $(round(data.enumeration_time, digits=2)) seconds")
    println("Created: $(data.timestamp)")
    
    # Cluster statistics
    total_clusters = sum(length(clusters) for clusters in values(data.clusters_by_site))
    println("\nCluster Statistics:")
    println("  Total clusters: $total_clusters")
    
    # Distribution by weight
    by_weight = Dict{Int, Int}()
    for clusters in values(data.clusters_by_site)
        for cluster in clusters
            by_weight[cluster.weight] = get(by_weight, cluster.weight, 0) + 1
        end
    end
    
    println("  Distribution by weight:")
    for weight in sort(collect(keys(by_weight)))
        count = by_weight[weight]
        println("    Weight $weight: $count clusters")
    end
    
    # Distribution by site
    cluster_counts = [length(clusters) for clusters in values(data.clusters_by_site)]
    println("  Distribution by site:")
    println("    Min: $(minimum(cluster_counts)) clusters")
    println("    Max: $(maximum(cluster_counts)) clusters") 
    println("    Mean: $(round(sum(cluster_counts) / length(cluster_counts), digits=2)) clusters")
    
    return nothing
end


## centralizing things
 

# Ursell function for connected clusters (from original implementation)
function ursell_function(cluster::Cluster, all_loops::Vector{Loop})
    """
    Compute the Ursell function φ(W) for a connected cluster W.
    For connected clusters, this is (-1)^(|W|-1) * (|W|-1)!
    where |W| is the total number of loops in the cluster (with multiplicities).
    """
    total_loops = cluster.total_loops
    if total_loops == 1
        return 1.0
    else
        return -0.5
    end
end

function cluster_contr(T_normalized, messages, edges, links, adj_mat, cluster, all_loops)
    # Compute Ursell function φ(W)
    phi_W = ursell_function(cluster, all_loops)
    
    # if abs(phi_W) < 1e-15
    #     return 0.0  # Skip negligible contributions
    # end
    
    # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
    Z_W = 1.0 + 0im  # Complex number for cluster contribution
    # invalidloop = false 
    for loop_id in cluster.loop_ids
        multiplicity = cluster.multiplicities[loop_id]
        loop = all_loops[loop_id]
        
        # Convert loop edges format for BP.jl (ensure v1 < v2 ordering)
        loop_edges_bp = Tuple{Int,Int}[]
        for edge in loop.edges
            v1, v2 = edge
            push!(loop_edges_bp, (min(v1, v2), max(v1, v2)))
        end
        
        # Compute loop contribution Z_l using BP
        try
            Z_l_tensor = BP.loop_contribution(loop_edges_bp, messages, T_normalized, edges, links, adj_mat)
            # Extract scalar value from ITensor
            Z_l = scalar(Z_l_tensor)
            Z_W *= Z_l^multiplicity
        catch e
            # println("⚠️  Error computing loop contribution for loop $loop_id: $e")
            # println("invalid loops present...", !(all([e in edges for e in loop_edges_bp])))
            return 0.0  # Return 0 instead of breaking
        end
    end
    
    # Contribution to log partition function (not free energy)
    # According to polymer theory: log(Z) = log(Z_BP) + Σ φ(W) Z_W
    contribution = phi_W * Z_W 
    return contribution
end 

function cluster_contr_by_site(cluster_data, TN_normalized, messages, edges, links, adj_mat)
    all_loops = cluster_data.all_loops  
    clusters_by_site = cluster_data.clusters_by_site

    clustercorrx = 0.0 + 0.0im  # Start with complex number

    # Deduplicate clusters manually (we know Set doesn't work)
    clusters_by_signature = Dict{Tuple, Vector{Tuple{Int, Cluster}}}()

    for site in 1:length(TN_normalized)
        clusters = clusters_by_site[site]
        for cluster in clusters
            signature = (cluster.weight, cluster.total_loops, sort(cluster.loop_ids), sort(collect(cluster.multiplicities)))
            if !haskey(clusters_by_signature, signature)
                clusters_by_signature[signature] = []
            end
            push!(clusters_by_signature[signature], (site, cluster))
        end
    end

    # Compute contributions for unique clusters
    for (sig, cluster_list) in clusters_by_signature
        cluster = cluster_list[1][2]  # Take first representative
        contribution = cluster_contr(TN_normalized, messages, edges, links, adj_mat, cluster, all_loops)
        if !isnan(contribution) && isfinite(contribution)
            clustercorrx += Complex(contribution)  # Ensure complex arithmetic
        end
    end
    return clustercorrx
end 

function load_latest_cluster_file(N, w;bc = "periodic", save_dir = "saved_clusters")
    if !isdir(save_dir)
        error("No saved_clusters directory found!")
    end
    files = readdir(save_dir)
    matching_files = filter(f -> endswith(f, ".jld2") && contains(f, "L$N") && contains(f, "w$w") && contains(f, "$bc"), files)
    if isempty(matching_files)
        error("No matching cluster file found for N=$N, w=$w in $save_dir!")
    end
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading cluster data from: $(latest_file)")
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end

# No module exports needed for include-style usage