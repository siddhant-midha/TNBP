using ITensors, ITensorMPS, Plots, LaTeXStrings
using ProgressMeter, Graphs, LinearAlgebra
using Combinatorics
using Statistics
push!(LOAD_PATH, "functions/")
using BP

function ortho(T::ITensor, v::Vector{Float64})
    indx = inds(T)
    n = length(indx)
    for idx in indx 
        vec = ITensor(v,idx)
        vec = vec / norm(vec)
        T = T - (T * vec) * vec 
        T = T / norm(T)
    end 
    return T 
end 

function controllable_tensor(i1_in, i2_in, i1_out, i2_out; η=0, orthog=true, noise=0, biasing=false)
    """
    Creates a controllable 4-index tensor with adjustable rank-1 structure and noise.
    
    This function generates a tensor that interpolates between a pure rank-1 tensor
    and a random tensor, allowing controlled studies of tensor network properties.
    
    Args:
        i1_in, i2_in, i1_out, i2_out: ITensor indices for the four legs of the tensor
        η: Mixing parameter controlling deviation from rank-1 structure (default: 0)
           - η=0: Pure rank-1 tensor  
           - η>0: Mix of rank-1 + random components
        orthog: Whether to orthogonalize random component against rank-1 part (default: true)
        noise: Amplitude of random noise added to the base vector (default: 0)
        biasing: If true, initializes with biased vector (ones + noise), else random (default: false)
    
    Returns:
        ITensor: A 4-index tensor with controllable structure
    """
    
    # Generate base vector for rank-1 component
    # biasing=true creates a vector favoring certain basis states
    vec = biasing ? (ones(dim(i1_in)) + noise * rand(dim(i1_in))) : rand(dim(i1_in))
    vec = vec / norm(vec)  # Normalize to unit length
    
    # Create rank-1 tensor: outer product of the same vector on all indices
    # This creates a fully separable tensor T[i,j,k,l] = v[i] * v[j] * v[k] * v[l]
    fp = ITensor(vec, i1_in) * ITensor(vec, i2_in) * ITensor(vec, i1_out) * ITensor(vec, i2_out)
    fp = fp / norm(fp)  # Normalize the rank-1 component
    
    # Create random component
    if orthog
        # Generate random tensor and orthogonalize it against the rank-1 vector
        # This ensures the random component is truly orthogonal to the rank-1 part
        fm = ortho(randomITensor(i1_in, i2_in, i1_out, i2_out), vec)
    else
        # Use completely random tensor (not orthogonalized)
        fm = randomITensor(i1_in, i2_in, i1_out, i2_out)
    end
    fm = fm / norm(fm)  # Normalize the random component
    
    # Linear combination: interpolate between rank-1 and random structures
    # η=0: pure rank-1 tensor
    # η→∞: dominated by random component
    f = fp + η * fm 
    
    # Note: Final normalization is commented out to preserve the mixing ratio
    f = f / norm(f)
    
    return f
end

function peps_controllable(N, T; η=0, ti=true, orthog=true, noise = 0, biasing = false)
    ## if ti = true all tensors are identical, else all random, different
    ## if return_peps = true returns the peps matrix, else returns the tensors list
    χ = 2
    vinds = [Index(χ, "v$(n)t$(t)") for n in 1:N, t in 1:T-1]
    hinds = [Index(χ, "n$(n)h$(t)") for n in 1:N-1, t in 1:T]

    down, up, left, right = Index(2,"down"), Index(2,"up"), Index(2,"left"), Index(2,"right")
    tens_main = controllable_tensor(down, up, left, right; η=η, orthog=orthog,noise=noise,biasing=biasing)

    tensors = []
    peps = Matrix{ITensor}(undef, T, N)
    for n = 1:N 
        for t = 1:T
            if ti 
                tens = copy(tens_main)
            else 
                down, up, left, right = Index(2,"down"), Index(2,"up"), Index(2,"left"), Index(2,"right")
                tens = controllable_tensor(down, up, left, right; η=η, orthog=orthog,noise=noise,biasing=biasing)
            end 
            # Vertical connections (up/down)
            if t == 1
                tens *= ITensor([1.,1.] ./ sqrt(2), down)
                tens *= delta(up, vinds[n, t])
            elseif t == T
                tens *= ITensor([1.,1.] ./ sqrt(2), up)
                tens *= delta(down, vinds[n, t - 1])
            else
                tens *= delta(down, vinds[n, t - 1])
                tens *= delta(up, vinds[n, t])
            end

            # Horizontal connections (left/right)
            if n == 1
                tens *= ITensor([1.,1.] ./ sqrt(2), left)
                tens *= delta(right, hinds[n, t])
            elseif n == N
                tens *= ITensor([1.,1.] ./ sqrt(2), right)
                tens *= delta(left, hinds[n - 1, t])
            else
                tens *= delta(left, hinds[n - 1, t])
                tens *= delta(right, hinds[n, t])
            end
            peps[t,n] = tens
            push!(tensors, tens)
        end
    end

    return tensors, peps
end


"""
    contract_peps_no_phys(peps::Matrix{ITensor}; cutoff=1E-8, maxdim=100)

Approximately contracts a 2D PEPS tensor network that has no physical indices.
The network is contracted to a single scalar value (its norm).

The PEPS is contracted row-by-row from top to bottom using the boundary MPS method.
The first row is contracted into an MPS. Then, each subsequent bulk row is
treated as an MPO and applied to the boundary MPS. Finally, the last row,
treated as an MPS, is contracted with the result.

Arguments:
- `peps`: A 2D array (Matrix) of ITensors representing the PEPS without physical legs.
- `cutoff`: The truncation error cutoff for SVD during MPO-MPS application.
- `maxdim`: The maximum allowed bond dimension for the boundary MPS.

Returns:
- A complex number representing the full contraction of the PEPS.
"""
function contract_peps_no_phys(peps::Matrix{ITensor}; cutoff=1E-8, maxdim=100)
    Ny, Nx = size(peps)
    if Ny < 2
        error("PEPS must have at least 2 rows for this contraction method.")
    end

    # --- 1. Contract the first row into an MPS ---
    # This MPS becomes the initial top boundary. Its site indices are the
    # downward-pointing virtual indices of the first row.
    # println("Contracting first row into boundary MPS...")
    boundary_mps = MPS(peps[1, :])

    # --- 2. Iteratively contract the bulk rows (2 to Ny-1) ---
    # println("Contracting bulk rows...")
    for i in 2:(Ny - 1)
        # println("  Applying MPO from row $i / $Ny")

        # The site indices of our current boundary_mps are the virtual indices
        # connecting the previously contracted part to the current row.
        top_links = siteinds(boundary_mps)

        # Construct an MPO from the PEPS row. We must prime the top virtual
        # links of the PEPS tensors so the MPO constructor correctly identifies
        # them as the "input" site indices.
        mpo_tensors = [prime(peps[i,j], top_links[j]) for j in 1:Nx]
        row_mpo = MPO(mpo_tensors)

        # Prime the boundary MPS to match the MPO's input site indices.
        boundary_mps_p = prime(boundary_mps)

        # Apply the MPO. The resulting MPS will have the MPO's output site
        # indices (the unprimed bottom_links of the PEPS row) as its new site indices.
        boundary_mps = apply(row_mpo, boundary_mps_p; cutoff=cutoff, maxdim=maxdim)
    end

    # --- 3. Contract with the final row ---
    # println("Contracting with final row...")
    # The last row of the PEPS is treated as an MPS.
    final_row_mps = MPS(peps[Ny, :])

    # The site indices of `boundary_mps` (from the bulk) and `final_row_mps`
    # are the same set of virtual indices connecting row Ny-1 and Ny.
    # The inner product gives the final scalar result.
    result = inner(boundary_mps, final_row_mps)

    # println("Contraction finished.")
    return (result)
end



function get_marginal(tensors,adj_mat,messages,index)
    nbrs = BP.get_nbrs(adj_mat, index)
    Z_local = tensors[index] 
    for nbr in nbrs
        Z_local *= messages[nbr,index] 
    end
    return Z_local
end 


function bp_contract(tensors; maxiter=200, annealing=0.9, normalise=true)
    """
    Performs belief propagation contraction on a tensor network.
    
    Args:
        tensors: Vector of ITensors representing the tensor network
        maxiter: Maximum number of BP iterations (default: 200)
        annealing: Damping factor for message updates (default: 0.9)
        normalise: Whether to normalize messages (default: true)
    
    Returns:
        Z: The contracted value (product of marginals)
    """
    
    adj_mat, edges, links = BP.get_adj_mat(tensors)
    messages = BP.get_messages(tensors, edges, links) 
    messages = BP.message_passing(tensors, messages, edges, links, adj_mat; 
                                 α=annealing, max_iters=maxiter, diagnose=false, normalise=normalise)
    
    marginals = [get_marginal(tensors, adj_mat, messages, index) for index = 1:length(tensors)]
    Z = prod(marginals)
    
    return scalar(Z)
end


# Parameters
N = 4
T = 4
η_range = 0.0001:0.01:0.5  # Test η from 0 to 1 in steps of 0.1
nsamples = 1000  # Number of samples per η value

# Fixed parameters
noise = 0
ti = true   
orthog = false   
normalise = true 

# Storage for results
η_vals = collect(η_range)
mean_errors = Float64[]
std_errors = Float64[]

println("Comparing contraction methods...")

for η in η_vals
    println("Testing η = $η")
    
    relative_errors = Float64[]
    
    for sample in 1:nsamples
        # Generate random PEPS
        tensors, peps = peps_controllable(N, T; η=η, ti=ti, orthog=orthog, noise=noise)
        
        # Exact contraction
        exact_result = contract_peps_no_phys(peps; cutoff=1E-9, maxdim=32)
        
        # BP contraction  
        bp_result = bp_contract(tensors; maxiter=500, annealing=0.9, normalise=normalise)
        
        # Compute relative error
        if abs(exact_result) > 1e-12  # Avoid division by very small numbers
            rel_error = abs(bp_result - exact_result) / abs(exact_result)
        else
            rel_error = abs(bp_result - exact_result)  # Absolute error if exact is ~0
        end
        
        push!(relative_errors, rel_error)
    end
    
    # Compute statistics
    push!(mean_errors, mean(relative_errors))
    push!(std_errors, std(relative_errors))
    
    println("  Mean relative error: $(mean_errors[end])")
    println("  Std relative error: $(std_errors[end])")
end

# Plot the results
p = plot(η_vals, mean_errors, 
         ribbon=std_errors,
         linewidth=2,
         grid=true,
         legend=:topleft,
         size=(600, 400))

# Add some styling
plot!(p)

display(p)

# Print summary statistics
println("\nSummary:")
println("η range: $(minimum(η_vals)) to $(maximum(η_vals))")
println("Samples per η: $nsamples")
println("Min mean error: $(minimum(mean_errors))")
println("Max mean error: $(maximum(mean_errors))")

readline() 