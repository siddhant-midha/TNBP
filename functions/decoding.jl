"""
Belief Propagation Tensor Network Decoding for Error Correction Codes

This module implements BP-based decoding using tensor networks for LDPC codes.
The tensor network structure follows:
- Tensors 1:n: Data bit tensors (variable nodes)  
- Tensors n+1:n+m: Parity check tensors (factor nodes)

The factor graph represents the decoding problem where data indices can be 
marginalized to compute bit-wise probability distributions.
"""

function tensorargmax(probs)
    """
    Extract the most probable value from a single-index probability tensor.
    
    Args:
        probs: ITensor with single index containing probability distribution
    
    Returns:
        Int: 0-indexed position of maximum probability (0, 1, ..., d-1)
    
    Note: Automatically normalizes probabilities and converts from Julia's 1-indexing.
    """
    ix = inds(probs)[1]
    probs = [real((probs)[ix=>n]) for n in 1:dim(ix)]
    probs ./= sum(probs)
    return Int(argmax(probs)-1) #sample_bit(probs[1])
end

function parity_tensor(index_arr, parity)
    """
    Create a parity constraint tensor for error correction codes.
    
    Args:
        index_arr: Array of ITensor indices representing variables in the constraint
        parity: Integer (0 or 1) specifying the required parity constraint
    
    Returns:
        ITensor: Constraint tensor with 1.0 for configurations satisfying the parity
                 and 0.0 elsewhere
    
    Note: Used in LDPC decoding where parity checks enforce sum(bits) % 2 == parity.
    """
    num = length(index_arr)
    tens = ITensor(index_arr)
    for i in 0:(2^num - 1)
        bits = digits(i, base=2, pad=num) 
        if sum(bits) % 2 == parity
            inds_tuple = (index_arr[j] => bits[j] + 1 for j in 1:num)
            tens[inds_tuple...] = 1.0        
        end 
    end
    return tens 
end 

function get_nbrs_of_check(pcmat, v)
    """
    Get data bit neighbors of a parity check node.
    
    Args:
        pcmat: Parity check matrix (rows=checks, columns=data bits)
        v: Check node index
    
    Returns:
        Vector{Int}: Indices of data bits connected to this check
    """
    row = pcmat[v, :]
    return findall(x -> x == 1, row)
end

function get_nbrs_of_data(pcmat, v)
    """
    Get parity check neighbors of a data bit node.
    
    Args:
        pcmat: Parity check matrix (rows=checks, columns=data bits)
        v: Data bit index
    
    Returns:
        Vector{Int}: Indices of checks connected to this data bit
    """
    col = pcmat[:, v]
    return findall(x -> x == 1, col)
end


function get_network(pcmat, syndrome, pbias)
    """
    Construct tensor network for LDPC decoding problem.
    
    Args:
        pcmat: Parity check matrix (m×n, where m=n-k)
        syndrome: Observed syndrome vector (length m)
        pbias: Bit flip probability for channel noise model
    
    Returns:
        Tuple of (data_tensors, syn_tensors, datainds):
        - data_tensors: Tensors for data bits with bias and connectivity
        - syn_tensors: Parity check constraint tensors 
        - datainds: Open indices for data bits (for marginal computation)
    
    Note: Creates a factor graph where data indices remain open for decoding.
    """
    ## pcmat: (n-k) x (n)
    ## syndrome: (n-k)
    ## close syndrome legs with given syndrome, keeps data legs open for later
    m, n = size(pcmat)
    k = n - m 
    indmat = [Index(2, "s$(i)d$(j)") for i in 1:m, j in 1:n]
    datainds = [Index(2, "x$i") for i in 1:n]
    data_tensors = []
    syn_tensors = []

    ## data data_tensors

    for i = 1:n 
        dummy = Index(2, "biasi")
        biastensor = ITensor(2 .* [1-pbias, pbias],dummy)
        checks = get_nbrs_of_data(pcmat,i)
        indxs = [indmat[jj,i] for jj in checks]
        push!(indxs,dummy)
        push!(indxs,datainds[i])
        push!(data_tensors,(delta(indxs) * biastensor))
    end 

    ## check tensors 
    for j = 1:m
        datas = get_nbrs_of_check(pcmat,j)
        tensor = parity_tensor([indmat[j,ii] for ii in datas],syndrome[j])
        push!(syn_tensors,tensor)
    end 
    return data_tensors, syn_tensors, datainds
end

function bit_to_onehot(x)
    if x == 0
        return [1, 0]
    elseif x == 1
        return [0, 1]
    elseif x == -1
        return [1, 1] ./ 2
    else
        error("Input must be -1, 0, or 1")
    end
end



function get_marginal_data_tensors(data_tensors,data_indices,data_inputs;exclude=[])
    """
    Marginalize data tensors by fixing observed bit values.
    
    Args:
        data_tensors: Vector of data bit tensors from get_network
        data_indices: Vector of open data indices 
        data_inputs: Vector of observed bit values (0, 1, or -1 for unknown)
        exclude: Indices to keep open (don't marginalize these bits)
    
    Returns:
        Vector of tensors with specified data legs closed/marginalized
    
    Note: Uses bit_to_onehot to convert observed bits to tensor form.
    """
    N = length(data_tensors)
    marginalized = []
    for i = 1 : N 
        if !(i in exclude)
            tens = data_tensors[i] *  ITensor(bit_to_onehot(data_inputs[i]),data_indices[i])
        else 
            tens = data_tensors[i]
        end 
        push!(marginalized,tens)
    end 
    return marginalized
end


function sample_bit(p0)
    """
    Sample a binary bit based on probability.
    
    Args:
        p0: Probability of returning 0
    
    Returns:
        Int: 0 with probability p0, 1 with probability (1-p0)
    """
    return rand() < p0 ? 0 : 1
end

function get_marginal(tensors,adj_mat,messages,index)
    """
    Compute marginal probability distribution for a tensor node.
    
    Args:
        tensors: Vector of ITensors in the network
        adj_mat: Adjacency matrix defining network connectivity
        messages: Message matrix from belief propagation
        index: Node index to compute marginal for
    
    Returns:
        ITensor: Marginal probability distribution for the specified node
    
    Note: Multiplies local tensor with normalized incoming messages from neighbors.
    """
    nbrs = BP.get_nbrs(adj_mat, index)
    Z_local = tensors[index] 
    for nbr in nbrs
        Z_local *= messages[nbr,index] / sqrt(scalar(messages[nbr,index] * messages[index,nbr])) 
    end
    return Z_local
end
