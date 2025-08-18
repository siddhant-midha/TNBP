using ITensors, ITensorMPS, Plots, LaTeXStrings
using ProgressMeter, Graphs, LinearAlgebra
using Combinatorics
using Statistics
using Test
using Serialization

include("dependencies.jl")
include("functions/ClusterEnumeration.jl")
include("functions/boundary_evolution.jl")

push!(LOAD_PATH, "functions/")
using BP

function load_latest_cluster_file(N,w)
    """Load the most recent cluster enumeration file."""
    save_dir = "saved_clusters"
    
    if !isdir(save_dir)
        error("No saved_clusters directory found!")
    end
    
    # Look for files matching our criteria (weight 10, size 11, PBC)
    files = readdir(save_dir)
    
    matching_files = filter(f -> contains(f, "L$N") && contains(f, "w$w"), files)
    
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    
    println("📖 Loading cluster data from: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    
    return loaded_data["data"]
end



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




function bp_contract(tensors, loop_dict; maxiter=500, annealing=0.9, normalise=true)
    """
    Performs belief propagation contraction on a tensor network with multiple loop orders.
    
    Args:
        tensors: Vector of ITensors representing the tensor network
        loop_dict: Dictionary with keys as loop orders and values as loop lists
                  (empty dict for standard BP only)
        maxiter: Maximum number of BP iterations (default: 500)
        annealing: Damping factor for message updates (default: 0.9)
        normalise: Whether to normalize messages (default: true)
    
    Returns:
        Dict with keys "vacuum" and loop orders, values are contracted results
    """
    adj_mat, edges, links = BP.get_adj_mat(tensors)
    messages = BP.get_messages(tensors, edges, links) 
    messages  = BP.message_passing(tensors, messages, edges, adj_mat; 
                                    α=annealing, max_iters=maxiter, diagnose=false, normalise=normalise)
    Z_list = BP.get_fixed_point_list(tensors,messages,adj_mat)    
    tensors = BP.normalize_tensors(tensors,Z_list)                      
    Z = prod(Z_list)
    
    # Initialize results dictionary
    results = Dict()
    results["vacuum"] = Z
    
    # If no loops provided, return only vacuum result
    if isempty(loop_dict)
        return results
    end
    
    # Compute contributions for each loop order
    for (order, loops) in loop_dict
        contribution = 0
        invalidloops = false 
        
        for loop in loops
            if all([e in edges for e in loop])
                contr = scalar(BP.loop_contribution(loop, messages, tensors, edges, links, adj_mat))
                contribution += contr
            else 
                invalidloops = true 
            end
        end 
        
        # if invalidloops
        #     println("Invalid loops were present for order $order, filtered out manually")
        # end 
        
        results[order] = Z * (1 + contribution)
    end
    
    return results
end

N = 10 
T = N
η_range = 0.0:0.2:2. 
nsamples = 100

w_list = [4,6,8,10]

# Load loop data for all weights
loop_data = Dict()
for w in w_list
    data = load_latest_cluster_file(N,w)
    loop_objects = data.all_loops
    loop_data[w] = [loop_object.edges for loop_object in loop_objects]
    println("Loaded $(length(loop_data[w])) loops for weight $w")
end

# Fixed parameters
noise = 0
ti = true   
orthog = false   
normalise = true 
annealing = .9 
maxiter = 500

# Storage for results - now including vacuum (no loops) and all weights
η_vals = collect(η_range)
mean_errors = Dict()
mean_errors["vacuum"] = Float64[]  # Standard BP without any loop corrections
for w in w_list
    mean_errors[w] = Float64[]
end

println("Comparing contraction methods across multiple loop orders...")

for η in η_vals
    println("Testing η = $η")
    
    # Storage for this η value
    relative_errors = Dict()
    relative_errors["vacuum"] = Float64[]
    for w in w_list
        relative_errors[w] = Float64[]
    end
    
    for sample in 1:nsamples
        # Generate random PEPS
        tensors, peps = peps_controllable(N, T; η=η, ti=ti, orthog=orthog, noise=noise)
        
        # Exact contraction
        exact_result = contract_peps_no_phys(peps; cutoff=1E-9, maxdim=32)
        
        # Single BP run with all loop corrections
        results = bp_contract(tensors, loop_data; maxiter=maxiter, annealing=annealing, normalise=normalise)
        
        # Compute relative errors for all methods
        for method in ["vacuum"; w_list]
            result = results[method]
            
            if abs(exact_result) > 1e-12
                rel_error = abs(result - exact_result) / abs(exact_result)
            else
                rel_error = abs(result - exact_result)
            end
            push!(relative_errors[method], rel_error)
        end
    end
    
    # Compute mean errors for this η
    push!(mean_errors["vacuum"], mean(relative_errors["vacuum"]))
    for w in w_list
        push!(mean_errors[w], mean(relative_errors[w]))
    end
    
    println("  Mean vacuum BP error: $(mean_errors["vacuum"][end])")
    for w in w_list
        println("  Mean w=$w loop-corrected error: $(mean_errors[w][end])")
    end
end

# Create comprehensive plot
p = plot(xlabel="η", 
         ylabel="Relative Error",
         title="BP vs Loop-Corrected PEPS Contraction Error (Multiple Orders)",
         yscale=:log10,
         grid=true,
         legend=:topleft,
         size=(800, 500))

# Define colors and markers for different methods
colors = [:blue, :red, :green, :orange, :purple]
markers = [:circle, :square, :diamond, :utriangle, :star5]

# Plot vacuum (standard BP)
plot!(p, η_vals, mean_errors["vacuum"],
      linewidth=2,
      markershape=markers[1],
      markersize=4,
      label="Vacuum (Standard BP)",
      color=colors[1])

# Plot each loop order
for (i, w) in enumerate(w_list)
    plot!(p, η_vals, mean_errors[w],
          linewidth=2,
          markershape=markers[i+1],
          markersize=4,
          label="Loop Order w=$w",
          color=colors[i+1])
end

# Add some styling
plot!(p, xlims=(minimum(η_vals)-0.01, maximum(η_vals)+0.01))

display(p)

# Print comprehensive summary statistics
println("\nSummary:")
println("η range: $(minimum(η_vals)) to $(maximum(η_vals))")
println("Samples per η: $nsamples")
println("\nVacuum (Standard BP):")
println("  Min mean error: $(minimum(mean_errors["vacuum"]))")
println("  Max mean error: $(maximum(mean_errors["vacuum"]))")

for w in w_list
    println("\nLoop Order w=$w:")
    println("  Number of loops: $(length(loop_data[w]))")
    println("  Min mean error: $(minimum(mean_errors[w]))")
    println("  Max mean error: $(maximum(mean_errors[w]))")
end

readline()