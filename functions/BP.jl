"""
BP.jl - Belief Propagation for Tensor Networks

This module implements belief propagation (BP) algorithms for computing partition functions
of tensor networks, with support for loop corrections via polymer/cluster expansion.

THEORY BACKGROUND:
-----------------
Belief propagation treats tensor network contraction as a message passing problem on
the dual graph, where:
- Each tensor becomes a vertex 
- Shared indices become edges
- Messages μ[i→j] approximate marginal distributions on edges
- Converged messages yield approximate partition function Z ≈ ∏ Z_local

For exact results, loop corrections account for correlations ignored by tree-like BP:
Z_exact = Z_BP * (1 + ∑ Loop_corrections + ∑ Higher_order_terms + ...)

MAIN WORKFLOW:
--------------
1. **Network Setup**: Extract graph structure from tensor network
   ```julia
   adj_mat, edges, links = get_adj_mat(tensors)
   ```

2. **Message Initialization**: Start with identity messages + noise
   ```julia
   messages = get_messages(tensors, edges, links; random_part=0.01)
   ```

3. **BP Iteration**: Converge messages via fixed-point iteration
   ```julia
   messages = message_passing(tensors, messages, edges, adj_mat; α=0.5, max_iters=1000)
   ```

4. **Partition Function**: Extract BP approximation
   ```julia
   Z_list = get_fixed_point_list(tensors, messages, adj_mat)
   Z_BP = prod(Z_list)
   ```

5. **Loop Corrections**: Add systematic improvements
   ```julia
   loop_correction = loop_contribution([(1,2), (2,3), (3,1)], messages, tensors, edges, links, adj_mat)
   Z_improved = Z_BP * (1 + scalar(loop_correction))
   ```

KEY FUNCTIONS:
--------------
**Network Analysis:**
- `get_adj_mat()`: Extract graph structure from tensor indices
- `get_nbrs()`: Find neighbors in adjacency matrix

**Message Passing:**  
- `get_messages()`: Initialize BP messages
- `message_passing()`: Iterative BP algorithm with damping/annealing
- `check_self_consistency()`: Validate message convergence

**Partition Functions:**
- `get_fixed_point_list()`: Local partition functions Z_v at each vertex
- `normalize_tensors()`: Normalize tensors by local Z values
- `mean_free_partition_fn()`: Partition function on vertex subsets

**Loop Corrections:**
- `excited_edge()`: Projector P⊥ = I - μ†μ orthogonal to messages  
- `loop_contribution()`: Systematic loop correction terms

TENSOR NETWORK CONVENTIONS:
---------------------------
- **Closed Networks**: No open/physical indices (for partition functions)
- **Index Sharing**: Connected tensors share exactly one index
- **Edge Ordering**: All edges stored as (i,j) with i < j
- **Message Direction**: μ[i,j] = message from vertex i to vertex j
- **Index Priming**: Used in loop corrections for proper contractions

ALGORITHMIC FEATURES:
--------------------
- **Damping/Annealing**: α < 1 for stable convergence in hard cases
- **Noise Injection**: Random perturbations prevent message decay
- **Normalization**: Messages kept at unit norm for numerical stability
- **Convergence Diagnostics**: Track message changes vs iterations
- **Self-Consistency Checks**: Validate BP fixed-point equations

PHYSICS APPLICATIONS:
--------------------
- **PEPS Contraction**: 2D tensor networks from quantum many-body systems
- **Statistical Mechanics**: Classical spin models, Ising systems
- **Quantum Circuits**: Amplitude estimation, sampling problems
- **Optimization**: QUBO, max-cut problems on graphs

PERFORMANCE NOTES:
------------------
- Complexity: O(χ³ × E × iterations) where χ = bond dimension, E = edges
- Memory: O(N × χ) for N tensors with bond dimension χ  
- Convergence: Typically 10-1000 iterations depending on problem difficulty
- Scaling: BP is approximate but scales polynomially vs exponential exact methods

EXAMPLE USAGE:
--------------
```julia
using ITensors
include("BP.jl")
using .BP

# Create 2×2 PEPS tensor network
N, T = 2, 2
tensors = create_peps_tensors(N, T)  # User-defined tensor creation

# Run BP analysis
adj_mat, edges, links = BP.get_adj_mat(tensors)
messages = BP.get_messages(tensors, edges, links)
messages = BP.message_passing(tensors, messages, edges, adj_mat; α=0.7, max_iters=500)

# Get results
Z_list = BP.get_fixed_point_list(tensors, messages, adj_mat)
Z_BP = prod(Z_list)
println("BP partition function: ", Z_BP)

# Add loop corrections
for edge_subset in combinations(edges, 1)  # Single-edge corrections
    correction = BP.loop_contribution(edge_subset, messages, tensors, edges, links, adj_mat)
    println("Loop correction: ", scalar(correction))
end
```

AUTHORS: [Your name/institution]
VERSION: 1.0
DEPENDENCIES: ITensors.jl
"""
module BP

using ITensors

# [Rest of the module code follows...]
function get_adj_mat(tensors)
    """
    Extract the adjacency matrix and edge information from a tensor network.
    
    This function analyzes the index structure of a tensor network to determine
    the connectivity graph and build the adjacency matrix needed for belief propagation.
    
    INPUTS:
    -------
    tensors : Vector{ITensor}
        List of ITensor objects representing the tensor network.
        Each tensor corresponds to a vertex in the graph.
        Tensors must have no open/external indices (closed tensor network).
        Connected tensors share exactly one common index.
    
    OUTPUTS:
    --------
    adj_mat : Matrix{Int32}
        Symmetric adjacency matrix of size n×n where n = length(tensors)
        adj_mat[i,j] = 1 if tensors i and j share a common index (are connected)
        adj_mat[i,j] = 0 if tensors i and j have no common indices
        Diagonal elements are always 0 (no self-loops)
        
    edges : Vector{Tuple{Int,Int}}
        List of all edges in the tensor network graph
        Each edge is represented as (i,j) where i < j (smaller vertex first)
        Length equals the number of shared indices in the tensor network
        
    links : Vector{Index}
        ITensor Index objects corresponding to each edge in the edges list
        links[k] is the shared index between tensors at edge edges[k]
        Used for message passing and tensor contractions
    
    ALGORITHM:
    ----------
    1. Initialize empty adjacency matrix and edge lists
    2. For each pair of tensors (i,j) with i < j:
       - Check if they share any common indices using commoninds()
       - If yes: mark as connected in adjacency matrix, add edge to list
       - Store the shared index for later use in BP message passing
    3. Enforce symmetry: adj_mat[i,j] = adj_mat[j,i] = 1 for connected pairs
    4. Assert exactly one shared index per connection (standard TN assumption)
    
    TENSOR NETWORK ASSUMPTIONS:
    ---------------------------
    - Each pair of connected tensors shares exactly one index
    - No open/dangling indices (closed tensor network for partition function)
    - Graph connectivity matches tensor index structure
    - Standard ITensor conventions for index naming and structure
    
    EXAMPLE:
    --------
    ```julia
    # 2×2 PEPS tensor network
    tensors = [T₁, T₂, T₃, T₄]  # 4 tensors in a square lattice
    adj_mat, edges, links = get_adj_mat(tensors)
    
    # Result:
    # adj_mat = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0]  # Square lattice connectivity
    # edges = [(1,2), (1,3), (2,4), (3,4)]              # Edges with i < j convention  
    # links = [idx₁₂, idx₁₃, idx₂₄, idx₃₄]               # Shared indices
    ```
    
    NOTES:
    ------
    - Critical for BP initialization: defines message passing structure
    - Edge ordering (i < j) convention used throughout BP module
    - Links array enables efficient tensor contractions during message updates
    - Adjacency matrix used for neighbor finding in BP iterations
    """
    ## CONVENTION: Each node in `tensors` is assigned a number 1:len(tensors), 
    ## and edges are stored (a,b) with a<b  
    n = length(tensors)
    adj_mat = zeros(Int32, n, n)  # Initialize n×n zero matrix
    edges = Tuple{Int,Int}[]      # Initialize edge list
    links = []                    # Initialize shared index list
    
    for i in 1:n
        for j in (i+1):n  # Only check upper triangle: avoid duplicates and self-loops
            shared_indices = commoninds(tensors[i], tensors[j])
            
            if !isempty(shared_indices)
                # Tensors i and j are connected - update adjacency matrix
                adj_mat[i,j] = 1
                adj_mat[j,i] = 1  # Symmetric connection
                
                # Add edge with standard ordering convention (i < j)
                push!(edges, (i,j))
                
                # Verify exactly one shared index (standard TN assumption)
                @assert length(shared_indices) == 1 "Expected exactly one common index between tensors $i and $j, found $(length(shared_indices))"
                
                # Store the shared index for message passing
                push!(links, shared_indices[1])
            end
        end
    end
    
    return adj_mat, edges, links
end


function get_nbrs(adj_mat, v)
    """
    Get the neighbors of vertex v from the adjacency matrix.
    
    Args:
        adj_mat: Adjacency matrix (n×n, symmetric, binary)
        v: Vertex index (1-based)
    
    Returns:
        Vector{Int}: List of neighbor vertex indices
    
    Example:
        adj_mat = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0]
        get_nbrs(adj_mat, 1) → [2, 3]
    """
    row = adj_mat[v, :]
    return findall(x -> x == 1, row)
end


function get_messages(tensors,edges,links;random_part=0)
    """
    Initialize BP messages for tensor network message passing.
    
    Args:
        tensors: Vector of ITensors in the network
        edges: List of (v1,v2) edges with v1 < v2 convention
        links: ITensor indices for each edge
        random_part: Amplitude of random noise added to identity messages (default: 0)
    
    Returns:
        Matrix{ITensor}: n×n message matrix μ[i,j] from vertex i to j
                        Initialized as normalized identity + noise on each shared index
    
    Notes:
        - Messages initialized as δ(index) + noise for numerical stability
        - Both directions μ[v1,v2] and μ[v2,v1] initialized independently
        - All messages normalized to unit norm
    """
    n = length(tensors)
    messages = Array{ITensor}(undef,n,n)
    for (e, edge) in enumerate(edges) 
        v1, v2 = edge 
        index = links[e]
        messages[v1, v2] = delta(ComplexF64, index) + random_part * randomITensor(index)
        messages[v2, v1] = delta(ComplexF64, index) + random_part * randomITensor(index)
        messages[v1, v2] = messages[v1, v2] / norm(messages[v1, v2])
        messages[v2, v1] = messages[v2, v1] / norm(messages[v2, v1])
    end 
    return messages
end


function check_self_consistency(tensors, messages, adj_mat; verbose=false)
    """
    Check the self-consistency of BP messages by verifying that each message
    equals the tensor contracted with all incoming messages except the one being checked.
    
    Args:
        tensors: Vector of ITensors in the tensor network
        messages: Matrix of converged BP messages
        adj_mat: Adjacency matrix of the tensor network
        verbose: If true, prints individual error contributions
    
    Returns:
        total_error: Sum of all message inconsistencies
        max_error: Maximum individual message error
    """
    
    total_error = 0.0
    max_error = 0.0
    
    for (v, T) in enumerate(tensors)
        nbrs = BP.get_nbrs(adj_mat, v)
        
        for n in nbrs 
            # Get the converged message from v to n
            converged = messages[v, n] 
            converged = converged / norm(converged)
            # Compute what this message should be: tensor v contracted with 
            # all incoming messages except the one from n
            selfconsistent = copy(T)
            for inp in setdiff(nbrs, n)
                selfconsistent *= messages[inp, v]
            end 
            selfconsistent = selfconsistent / norm(selfconsistent)
            # Compute the error for this message
            err_vn = norm(converged - selfconsistent)
            
            if verbose
                println("Message ($v → $n): error = $err_vn")
            end
            
            total_error += err_vn
            max_error = max(max_error, err_vn)
        end 
    end 
    
    return total_error, max_error
end


function message_passing(tensors,messages,edges,adj_mat;α=1,noise=0,max_iters=1000,diagnose=false,normalise=true)
    """
    Iterative belief propagation message passing for tensor networks.
    
    Args:
        tensors: Vector of ITensors in the network
        messages: Initial message matrix from get_messages()
        edges: List of (v1,v2) edges with v1 < v2 convention
        adj_mat: Adjacency matrix of the tensor network
        α: Damping factor (default: 1). α=1 instant updates, α<1 annealed
        noise: Random noise amplitude added each iteration (default: 0)
        max_iters: Maximum iterations before stopping (default: 1000)
        diagnose: If true, returns convergence history (default: false)
        normalise: Whether to normalize messages each iteration (default: true)
    
    Returns:
        Matrix{ITensor}: Converged messages μ[i,j] from vertex i to j
        Vector{Float64}: Convergence history (only if diagnose=true)
    
    Notes:
        - Updates messages bidirectionally: v1→v2 then v2→v1 for each edge
        - Converges when ||Δmessages|| < 1e-12 or max_iters reached
        - Message update: μ[v→w] = T[v] * ∏[neighbors u≠w] μ[u→v]
    """
    Δ = 100 
    arr = []
    iters = 0
    while Δ > 1e-12 && iters < max_iters
        iters += 1 
        δ = 0 
        
        # Forward pass: update v1 → v2 for each edge
        for e in edges 
            v1, v2 = e 
            @assert norm(messages[v1,v2]) > 1e-12 "messages decaying to zero!"
            update = tensors[v1]
            for nbr in get_nbrs(adj_mat, v1)
                if nbr != v2 
                    update = update * messages[nbr,v1]
                end
            end 
            new_message = normalise ? update/norm(update) : update 
            δ += norm(messages[v1,v2] - new_message)
            inde = inds(messages[v1,v2])[1]
            messages[v1,v2] = (1-α) * messages[v1,v2] + α * new_message 
                                    + noise * ITensor(randn(dim(inde)),inde)
            messages[v1,v2] = normalise ? messages[v1,v2] / norm(messages[v1,v2]) : messages[v1,v2]
        end 
        
        # Backward pass: update v2 → v1 for each edge
        for e in edges
            v2, v1 = e 
            @assert norm(messages[v1,v2]) > 1e-12 "messages decaying to zero!"
            update = tensors[v1]
            for nbr in get_nbrs(adj_mat, v1)
                if nbr != v2 
                    update = update * messages[nbr,v1] 
                end
            end 
            new_message = normalise ? update/norm(update) : update 
            δ += norm(messages[v1,v2] - new_message)
            inde = inds(messages[v1,v2])[1]
            messages[v1,v2] = (1-α) * messages[v1,v2] + α * new_message 
                                    + noise * ITensor(randn(dim(inde)),inde)
            messages[v1,v2] = normalise ? messages[v1,v2] / norm(messages[v1,v2]) : messages[v1,v2]
        end   
        
        Δ = δ  
        push!(arr,δ)
    end 
    
    if diagnose 
        return messages, arr
    else 
        return messages 
    end 
end

function excited_edge(edge,messages,edges,links)
    """
    Compute the excited projector P⊥ = I - μ†μ for a tensor network edge.
    
    Args:
        edge: Edge tuple (v1,v2) to excite
        messages: BP message matrix μ[i,j] from vertex i to j
        edges: List of all edges with v1 < v2 convention
        links: ITensor indices for each edge
    
    Returns:
        ITensor: Excited projector P⊥ = I - μ†_{v1→v2} * μ_{v2→v1}
    
    Notes:
        - Projects onto subspace orthogonal to BP messages
        - Used in loop corrections: replaces standard edge with excited state
        - Primes v1→v2 message to match tensor contraction conventions
    """
    v1, v2 = edge
    v1, v2 = min(v1,v2), max(v1,v2)
    edge_index = findfirst(e -> e == edge, edges)
    index = links[edge_index]
    iden = ITensor(index,prime(index)) 
    for n in 1:dim(index)
        iden[index=>n, prime(index)=>n] = 1.0 + 0im  # Fill diagonal with ComplexF64 1s
    end
    return iden - (prime(messages[v1,v2]) * messages[v2,v1]) / (messages[v1,v2] * messages[v2,v1]) 
                                            ## convention: in the excited edge, the daggered index is v1 -> v2 with v1 < v2
                                            ## so the corresponding index on T[v2] needs to be primed
end




function loop_contribution(loop, messages, tensors, edges, links, adj_mat)
    """
    Compute the loop correction tensor for the tensor network partition function.
    
    This implements the mathematical formula from polymer/cluster expansion:
    Z_l = (∏[loop edges] P⊥_{edge}) * (∏[non-loop edges] μ_{edge}) * (∏[all vertices] T_v)
    
    where P⊥_{edge} = I - μ†_{v1→v2} * μ_{v2→v1} is the excited projector.
    
    INPUTS:
    -------
    loop : Vector{Tuple{Int,Int}}
        List of edges forming the loop, e.g., [(1,2), (2,3), (3,1)]
        CRITICAL: Each edge (v1,v2) must satisfy v1 < v2 (smaller vertex first)
        
    messages : Matrix{ITensor}
        BP messages μ[i,j] from vertex i to vertex j
        Size: n×n where n = number of vertices
        
    tensors : Vector{ITensor}
        Original tensor network tensors T_v at each vertex
        
    edges : Vector{Tuple{Int,Int}}
        Complete list of all edges in tensor network, with v1 < v2 convention
        
    links : Vector{Index}
        ITensor indices corresponding to each edge
        links[i] is the shared index for edge edges[i]
        
    adj_mat : Matrix{Int32}
        Adjacency matrix of the tensor network graph
        adj_mat[i,j] = 1 if vertices i,j are connected, 0 otherwise
    
    OUTPUT:
    -------
    ITensor : The loop correction tensor
        If all indices contract to scalar: use scalar() to extract ComplexF64 value
        If open indices remain: tensor can be further contracted with other objects
        
    ALGORITHM:
    ----------
    1. Replace loop edge messages with excited projectors P⊥ = I - μ†μ
    2. For each vertex in loop:
       - Apply appropriate index priming to match projector structure  
       - Contract vertex tensor T_v with projectors
       - Add BP messages from non-loop neighbors
    3. Return contracted tensor (may be scalar or have remaining open indices)
    
    PHYSICS:
    --------
    This computes how the partition function changes when "exciting" a loop
    - Loop edges use projectors orthogonal to BP messages
    - Non-loop edges use standard BP messages  
    - Result gives the loop correction term needed for tensor network loop expansion
    
    USAGE:
    ------
    ```julia
    # For closed loops (no open indices)
    loop_tensor = loop_contribution([(1,2), (2,3), (3,1)], messages, tensors, edges, links, adj_mat)
    loop_scalar = scalar(loop_tensor)  # Extract scalar value
    
    # For partial contractions (may have open indices)
    loop_tensor = loop_contribution([(1,2)], messages, tensors, edges, links, adj_mat)
    # Contract further or check inds(loop_tensor) for remaining indices
    ```
    """
    # Initialize tracking variables
    vertices_done = Set()  # Track which vertices have been processed
    loop_contri = 1        # Accumulate the loop contribution tensor product
    N = length(tensors)    # Total number of vertices in tensor network
    
    # Step 1: Process each edge in the loop
    for edge in loop
        v1, v2 = edge  # v1 < v2 by convention 
        
        # Replace BP message with excited projector P⊥ = I - μ†μ
        excitation = excited_edge(edge,messages,edges,links)
        loop_contri *= excitation
        
        # Step 2a: Process vertex v1 (if not already done)
        if !(v1 in vertices_done)
            vertices_done = union(vertices_done, v1)
            
            # Find all vertices connected to v1 within the loop
            excited_neighbors = Set([other_vertex for edge in loop if v1 in edge for other_vertex in edge if other_vertex != v1])
            
            # Determine which indices need priming for proper tensor contraction
            edges_with_v1 = filter(t -> v1 in t, loop)  # Find loop edges containing v1
            edge_indices = [findfirst(isequal(t), edges) for t in edges_with_v1]  # Map to global edge list
            larger_bools = [v1 == max(t[1], t[2]) for t in edges_with_v1]  # v1 is larger vertex → prime needed
            selected_links = links[edge_indices[larger_bools]]  # Select indices to prime
            
            # Apply tensor T_v1 with appropriate index priming
            contri = !isempty(selected_links) ? prime(tensors[v1], selected_links...) : tensors[v1]
            loop_contri *= contri
            
            # Add BP messages from neighbors outside the loop
            for w in get_nbrs(adj_mat, v1)
                current_edge = (min(v1,w), max(v1,w))  # Standard edge ordering
                edge_index = findfirst(e -> e == current_edge, edges)
                current_link = links[edge_index]
                
                # Only add message if: (1) index exists in current contraction, (2) neighbor is outside loop
                if current_link in inds(loop_contri) && !(w in excited_neighbors)
                    loop_contri *= messages[w,v1]  
                end
            end
        end
        
        # Step 2b: Process vertex v2 (if not already done) - identical logic to v1
        if !(v2 in vertices_done)
            vertices_done = union(vertices_done, v2)
            excited_neighbors = Set([other_vertex for edge in loop if v2 in edge for other_vertex in edge if other_vertex != v2])
            edges_with_v2 = filter(t -> v2 in t, loop) 
            edge_indices = [findfirst(isequal(t), edges) for t in edges_with_v2]
            larger_bools = [v2 == max(t[1], t[2]) for t in edges_with_v2]
            selected_links = links[edge_indices[larger_bools]]
            contri = !isempty(selected_links) ? prime(tensors[v2], selected_links...) : tensors[v2]
            loop_contri *= contri 

            for w in get_nbrs(adj_mat, v2)
                current_edge = (min(v2,w), max(v2,w))  
                edge_index = findfirst(e -> e == current_edge, edges)
                current_link = links[edge_index]

                if current_link in inds(loop_contri) && !(w in excited_neighbors)
                    loop_contri *= messages[w,v2]  
                end
            end
        end
    end    
    return loop_contri
end

function get_fixed_point_list(tensors,messages,adj_mat)
    """
    Compute local BP partition functions at each vertex.
    
    Args:
        tensors: Vector of ITensors in the network
        messages: Converged BP message matrix μ[i,j]
        adj_mat: Adjacency matrix of the tensor network
    
    Returns:
        Vector{ComplexF64}: Local partition function Z_v at each vertex
    
    Notes:
        - Z_v = T_v * ∏[neighbors] μ[nbr→v] / √(μ[nbr→v] * μ[v→nbr])
        - Used for tensor normalization and global partition function
    """
    Z_list = []
    for index = 1:length(tensors)
        nbrs = get_nbrs(adj_mat, index)
        Z_local = tensors[index] 
        for nbr in nbrs
            Z_local *= messages[nbr,index] / sqrt(scalar(messages[nbr,index] * messages[index,nbr])) 
        end
        @assert isempty(inds(Z_local))  "T[$index] must be a scalar"
        push!(Z_list,scalar(Z_local))
    end
    return Z_list  
end

function normalize_tensors(tensors,Z_list)
    """
    Normalize tensor network by local BP partition functions.
    
    Args:
        tensors: Vector of ITensors in the network
        Z_list: Local partition functions from get_fixed_point_list()
    
    Returns:
        Vector{ITensor}: Normalized tensors T'_v = T_v / Z_v
    
    Notes:
        - Creates normalized tensors for improved numerical stability
        - Preserves tensor network structure and connectivity
    """
    T_copy = copy(tensors)
    for index = 1:length(tensors)
        T_copy[index] /= Z_list[index]
    end
    return T_copy
end

function mean_free_partition_fn(these_vertices,tensors,messages,adj_mat)
    """
    Compute BP partition function on subset of vertices.
    
    Args:
        these_vertices: Set/collection of vertex indices to include
        tensors: Vector of ITensors in the network
        messages: Converged BP message matrix μ[i,j]
        adj_mat: Adjacency matrix of the tensor network
    
    Returns:
        ComplexF64: Product of local partition functions Z = ∏[v ∈ subset] Z_v
    
    Notes:
        - these_vertices = 1:N gives full BP partition function
        - Used in loop corrections for partial contractions
    """
    Z_list = get_fixed_point_list(tensors,messages,adj_mat)
    return prod([Z_list[v] for v in these_vertices])
end

end