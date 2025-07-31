module BP 

using ITensors

#= 
Expectation: A list `tensors` comprising the tensors in the TN is provided, with no open indices, and correctly mapped out common indices between different tensors in the list. 
=#

function get_adj_mat(tensors)
    #= get adjacency matrix, list of edges, 
    list of indices corresponding to edges of the TN `tensors` =# 
    ## CONVENTION: Each node in `tensors` is assigned a number 1:len(tensors), 
    ## and edges are stored (a,b) with a<b  
    n = length(tensors)
    adj_mat = zeros(Int32, n, n)  # Initialize n x n zero matrix
    edges = Tuple{Int,Int}[]      # Initialize edge list
    links = []
    for i in 1:n
        for j in (i+1):n  # Avoid duplicate checks and self-loops
            if !isempty(commoninds(tensors[i], tensors[j]))
                adj_mat[i,j] = 1
                adj_mat[j,i] = 1  # Symmetric connection
                push!(edges, (i,j)) ## i < j convention
                @assert length(commoninds(tensors[i], tensors[j])) == 1 "Expected to have only one common index."
                push!(links,commoninds(tensors[i], tensors[j])[1])
            end
        end
    end
    return adj_mat, edges, links
end


function get_nbrs(adj_mat, v)
    ## get nbrs of node v from adj_mat
    row = adj_mat[v, :]
    return findall(x -> x == 1, row)
end


function get_messages(tensors,edges,links;random_part=0)
    ## initialize the messages, 
    ## tensors: TN as a list 
    ## links: list of indices corresponding to edges in the tensor network
    ## edges: list of  (v1,v2) edges with convention v1 < v2 
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


function message_passing(tensors,messages,edges,links,adj_mat;α=1,noise=0,max_iters=1000,diagnose=false,normalise=true)
    ## tensors: tensor network 
    ## messages: matrix of message tensors {μ[v1,v2]}
    ## links: list of indices corresponding to edges in the tensor network
    ## edges: list of (v1,v2) edges with convention v1 < v2 
    ## α: damping controls message update rate, α=1 updates instantly, α < 1 anneals slowly 
    Δ = 100 
    arr = []
    iters = 0
    while Δ > 1e-3 && iters < max_iters
        iters += 1 
        δ = 0 
        for e in edges 
            
            ## want to update the message from v1 → v2 
            v1, v2 = e 
            @assert norm(messages[v1,v2]) > 1e-12 "messages decaying to zero!"
            update = tensors[v1]
            for nbr in get_nbrs(adj_mat, v1)
                if nbr != v2 
                    update = update * messages[nbr,v1]
                end
            end 
            # println(check_permutation_indices(messages[v1,v2],update))
            new_message = normalise ? update/norm(update) : update 
            δ += norm(messages[v1,v2] - new_message)
            inde = inds(messages[v1,v2])[1]
            messages[v1,v2] = (1-α) * messages[v1,v2] + α * new_message 
                                    + noise * ITensor(randn(dim(inde)),inde)
            messages[v1,v2] = normalise ? messages[v1,v2] / norm(messages[v1,v2]) : messages[v1,v2]

        end 
        
        # backward
        
        for e in edges
            
            ## want to update the message from v1 → v2 
            v2, v1 = e 
            @assert norm(messages[v1,v2]) > 1e-12 "messages decaying to zero!"
            update = tensors[v1]
            for nbr in get_nbrs(adj_mat, v1)
                if nbr != v2 
                    update = update * messages[nbr,v1] 
                end
            end 
            # println(check_permutation_indices(messages[v1,v2],update))
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


function mean_free_partition_fn(these_vertices,tensors,messages,edges,links,adj_mat)
    ## computes the fixed point partition function at the vertices specified by these_vertices 
    ## other things as usual, these_vertices = 1:N gives the BP fixed point
    Z = 1 
    for index = these_vertices
        nbrs = get_nbrs(adj_mat, index)
        Z_local = tensors[index] 
        for nbr in nbrs
            Z_local *= messages[nbr,index] 
        end
        @assert isempty(inds(Z_local))  "tensors[$index] must be a scalar"
        Z *= scalar(Z_local)
    end
    return Z # -log(Z)/(L^2)
end 


function excited_edge(edge,messages,edges,links,adj_mat)
    ## gives the excited projector on the specified edge
    v1, v2 = edge
    v1, v2 = min(v1,v2), max(v1,v2)
    edge_index = findfirst(e -> e == edge, edges)
    index = links[edge_index]
    iden = ITensor(index,prime(index)) 
    for n in 1:dim(index)
        iden[index=>n, prime(index)=>n] = 1.0 + 0im  # Fill diagonal with ComplexF64 1s
    end
    return iden - prime(messages[v1,v2]) * messages[v2,v1] ## convention: in the excited edge, the daggered index is v1 -> v2 with v1 < v2
                                                           ## so the corresponding index on T[v2] needs to be primed
end 


function loop_contribution(loop, messages, tensors, edges, links, adj_mat)
    """
    Compute the loop correction Z_l to the tensor network partition function.
    
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
    ComplexF64 : The loop correction Z_l
        Complex scalar representing the contribution of this loop to log(Z)
        
    ALGORITHM:
    ----------
    1. Replace loop edge messages with excited projectors P⊥ = I - μ†μ
    2. For each vertex in loop:
       - Apply appropriate index priming to match projector structure  
       - Contract vertex tensor T_v with projectors
       - Add BP messages from non-loop neighbors
    3. Contract remaining vertices outside loop using standard BP
    4. Return scalar result
    
    PHYSICS:
    --------
    This computes how the partition function changes when "exciting" a loop
    - Loop edges use projectors orthogonal to BP messages
    - Non-loop edges use standard BP messages  
    - Result gives the loop correction term needed for tensor network loop expansion
    """
    # Initialize tracking variables
    vertices_done = Set()  # Track which vertices have been processed
    loop_contri = 1        # Accumulate the loop contribution tensor product
    N = length(tensors)    # Total number of vertices in tensor network
    
    # Step 1: Process each edge in the loop
    for edge in loop
        v1, v2 = edge  # v1 < v2 by convention 
        
        # Replace BP message with excited projector P⊥ = I - μ†μ
        excitation = excited_edge(edge, messages, edges, links, adj_mat)
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
    
    # Step 3: Contract to scalar and multiply by BP partition function on remaining vertices
    # Loop contribution (scalar) × BP partition function on vertices outside the loop
    return scalar(loop_contri) * mean_free_partition_fn(setdiff(Set(1:N), vertices_done), tensors, messages, edges, links, adj_mat)
end

function get_fixed_point_list(tensors,messages,edges,links,adj_mat)
    ## get a list of the fixed point partition function at each vertex
    Z_l = []
    for index = 1:length(tensors)
        nbrs = BP.get_nbrs(adj_mat, index)
        Z_local = tensors[index] 
        for nbr in nbrs
            Z_local *= messages[nbr,index] 
        end
        @assert isempty(inds(Z_local))  "T[$index] must be a scalar"
        push!(Z_l,scalar(Z_local))
    end
    return Z_l  
end

function normalize_tensor(tensors,Z_l)
    ## normalize the tensor network by dividing by the BP fixed point
    T_copy = copy(tensors)
    for index = 1:length(tensors)
        T_copy[index] /= Z_l[index]
    end
    return T_copy
end

end