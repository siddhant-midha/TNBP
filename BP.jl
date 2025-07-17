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


function loop_contribution(loop,messages,tensors,edges,links,adj_mat)
    ## gives the excited loop contribution to the partition function at the specified loop, and BP solution at all other edges
    ## loop looks like [(v1,v2),(v2,v3) ...] ## MAKE SURE all edges in the loop are oriented such that v1 < v2
    vertices_done = Set()
    loop_contri = 1 
    N = length(tensors)
    for edge in loop
        v1, v2 = edge ### v1 < v2 by convention 
        excitation = excited_edge(edge,messages,edges,links,adj_mat)
        loop_contri *= excitation
        if !(v1 in vertices_done)
            vertices_done = union(vertices_done, v1)
            excited_neighbors = Set([other_vertex for edge in loop if v1 in edge for other_vertex in edge if other_vertex != v1])
            edges_with_v1 = filter(t -> v1 in t, loop)  ## find excited edges with v1
            edge_indices = [findfirst(isequal(t), edges) for t in edges_with_v1] ## find the corresponding indices in edges list 
            larger_bools = [v1 == max(t[1], t[2]) for t in edges_with_v1] ## check if v1 is the larger vertex in the edge, so that we can prime it 
            selected_links = links[edge_indices[larger_bools]] ## thus now we select the links to be primed 
            contri = !isempty(selected_links) ? prime(tensors[v1], selected_links...) : tensors[v1] ## prime the selected links
            loop_contri *= contri
            for w in get_nbrs(adj_mat, v1)
                current_edge = (min(v1,w), max(v1,w))  
                edge_index = findfirst(e -> e == current_edge, edges)
                current_link = links[edge_index]
                if current_link in inds(loop_contri) && !(w in excited_neighbors)
                    loop_contri *= messages[w,v1]  
                end
            end
        end
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
    return scalar(loop_contri) * mean_free_partition_fn(setdiff(Set(1:N),(vertices_done)),tensors,messages,edges,links,adj_mat)
end

end 