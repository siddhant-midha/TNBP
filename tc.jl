push!(LOAD_PATH, "functions/")
using BP
using Random, Plots, SparseArrays, ITensors, Statistics, ProgressMeter, Colors, LinearAlgebra, CSV, DataFrames, Dates
include("ldpc_tanner_loops.jl")
include("functions/tc_decode.jl")


function parity_tensor(index_arr, parity)
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

function get_nbrs_of_check(adj_mat, v)
    ## get nbrs of node v from adj_mat
    row = adj_mat[v, :]
    return findall(x -> x == 1, row)
end

function get_nbrs_of_data(adj_mat, v)
    ## get nbrs of node v from adj_mat
    col = adj_mat[:, v]
    return findall(x -> x == 1, col)
end

function get_network(pcmat, syndrome, pbias)
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
        return [2, 0]
    elseif x == 1
        return [0, 2]
    elseif x == -1
        return [1, 1]
    else
        error("Input must be -1, 0, or 1")
    end
end


function get_marginal_data_tensors(data_tensors,data_indices,data_inputs;exclude=[])
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
    return rand() < p0 ? 0 : 1
end


function toric_code_X_parity_matrix(L::Int)
    N = 2 * L^2   # number of qubits
    M = L^2       # number of X checks (vertices)
    
    pcmat = spzeros(Int, M, N)
    
    for v in 1:M
        i = div(v-1, L) + 1  # vertex row
        j = mod1(v, L)       # vertex col
        
        # Horizontal edges (apply periodic boundary for i-1)
        h_edge1 = (mod1(i-1, L) - 1)*L + j
        h_edge2 = (mod1(i, L) - 1)*L + j
        
        # Vertical edges
        v_edge1 = L^2 + (i-1)*L + mod1(j-1, L)
        v_edge2 = L^2 + (i-1)*L + j
        
        pcmat[v, h_edge1] = 1
        pcmat[v, h_edge2] = 1
        pcmat[v, v_edge1] = 1
        pcmat[v, v_edge2] = 1
    end
    
    return pcmat
end

function get_marginal(tensors,adj_mat,messages,index)
    nbrs = BP.get_nbrs(adj_mat, index)
    Z_local = tensors[index] 
    for nbr in nbrs
        Z_local *= messages[nbr,index] 
    end
    return Z_local
end 

function loop_contribution(loop, messages, tensors, edges, links, adj_mat)
    # Initialize tracking variables
    vertices_done = Set()  # Track which vertices have been processed
    loop_contri = 1        # Accumulate the loop contribution tensor product
    N = length(tensors)    # Total number of vertices in tensor network
    
    # Step 1: Process each edge in the loop
    for edge in loop
        v1, v2 = edge  # v1 < v2 by convention 
        
        # Replace BP message with excited projector P⊥ = I - μ†μ
        excitation = BP.excited_edge(edge, messages, edges, links, adj_mat)
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
            for w in BP.get_nbrs(adj_mat, v1)
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

            for w in BP.get_nbrs(adj_mat, v2)
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
    return (loop_contri) #* mean_free_partition_fn(setdiff(Set(1:N), vertices_done), tensors, messages, edges, links, adj_mat)
end

function decode(pcmat, p, numsamples, L; pbias = 0.1, max_loop_order = 8)
    ## keeps a list of priors on all bits, and adds to the prior as and when a loop is encountered
    ## should work for arbitrary loop orders
    m, n = size(pcmat)
    tannerloopslist = [find_tanner_loops(pcmat, d; max_length=max_loop_order) for d in 1:n]
    # Track errors within this batch for both methods
    logical_errors_loops = zeros(numsamples)
    failure_rate_loops = zeros(numsamples)
    logical_errors_no_loops = zeros(numsamples)
    failure_rate_no_loops = zeros(numsamples)
    
    # Create progress bar for samples
    sample_prog = Progress(numsamples, desc="Samples: ", showspeed=true)
    
    for samp in 1:numsamples
        # Sample errors iid with probability p
        errors_true = [sample_bit(1-p) for _ in 1:n]
        # Compute syndrome
        syndrome = pcmat * errors_true .% 2
        # Decode
        data_tensors, syn_tensors, data_indices = get_network(pcmat, syndrome, pbias)
        errors_loops = Int.(-1 .* ones(n))
        errors_no_loops = Int.(-1 .* ones(n))
        tensors = vcat(get_marginal_data_tensors(data_tensors, data_indices, errors_loops), syn_tensors)
        adj_mat, edges, links = BP.get_adj_mat(tensors)
        messages = BP.get_messages(tensors,edges,links) 
        messages = BP.message_passing(tensors,messages,edges,links,adj_mat;α=0.95, max_iters=500,diagnose=false,normalise=true)
        
        # Two sets of priors: with and without loop corrections
        priors_loops = [ITensor([0.,0.], data_indices[i]) for i in 1:n]
        priors_no_loops = [ITensor([0.,0.], data_indices[i]) for i in 1:n]

        for d = 1:n
            ## vacuum contribution (same for both methods)
            probs = get_marginal(vcat(data_tensors,syn_tensors),adj_mat,messages,d)
            priors_loops[d] += probs
            priors_no_loops[d] += probs
            
            # Decode with loop corrections
            ix = inds(priors_loops[d])[1]
            probs_loops = [real((priors_loops[d])[ix=>n]) for n in 1:dim(ix)]
            probs_loops ./= sum(probs_loops)
            error_i_loops = sample_bit(probs_loops[1])
            errors_loops[d] = error_i_loops
            
            # Decode without loop corrections (use same base probabilities)
            ix_no = inds(priors_no_loops[d])[1]
            probs_no_loops = [real((priors_no_loops[d])[ix_no=>n]) for n in 1:dim(ix_no)]
            probs_no_loops ./= sum(probs_no_loops)
            error_i_no_loops = sample_bit(probs_no_loops[1])
            errors_no_loops[d] = error_i_no_loops

            # Early termination check for loop method
            if errors_loops[d] != errors_true[d]
                # no point in going forward for loop method
                break
            end 

            # Loop corrections (only for the loop method)
            tannerloops = tannerloopslist[d]
            loop_list = [tannerloop.edges for tannerloop in tannerloops] 
            data_bits_involved_list = [tannerloop.data_bits for tannerloop in tannerloops]
            check_bits_involved_list = [tannerloop.check_bits for tannerloop in tannerloops]
            
            for (i, loop) in enumerate(loop_list)
                data_bits_involved = data_bits_involved_list[i] 
                check_bits_involved = check_bits_involved_list[i] 
                data_bits_involved_other = collect(setdiff(data_bits_involved, [d]) )
                
                for data_bit in data_bits_involved_other
                    if true #all([errors_loops[bit] != -1 for bit in collect(setdiff(data_bits_involved, [data_bit]))])
                        mtensors = vcat(get_marginal_data_tensors(data_tensors, data_indices, errors_loops; exclude=[data_bit]), syn_tensors)
                        ## account for normalization from other nodes in the loop 
                        normlz = (prod([get_marginal(mtensors,adj_mat,messages,other_data_bit) for other_data_bit in collect(setdiff(data_bits_involved, [data_bit]) )]))
                        normlz *= (prod([get_marginal(mtensors,adj_mat,messages,n+check_bit) for check_bit in check_bits_involved]))
                        ## update the prior of the data bit in the loop (only for loop method)
                        priors_loops[data_bit] += loop_contribution(loop, messages, mtensors, edges, links, adj_mat) / normlz
                    end 
                end 
            end 
        end
        
        # Check results for both methods
        (syndrome_zero_loops, homology_trivial_loops) = check_decode_edges(errors_loops, errors_true, L)
        logical_errors_loops[samp] = (syndrome_zero_loops && homology_trivial_loops) ? 0.0 : 1.0
        failure_rate_loops[samp] = syndrome_zero_loops ? 0.0 : 1.0
        
        (syndrome_zero_no_loops, homology_trivial_no_loops) = check_decode_edges(errors_no_loops, errors_true, L)
        logical_errors_no_loops[samp] = (syndrome_zero_no_loops && homology_trivial_no_loops) ? 0.0 : 1.0
        failure_rate_no_loops[samp] = syndrome_zero_no_loops ? 0.0 : 1.0
        
        # Update sample progress
        next!(sample_prog)
    end
    
    return (mean(logical_errors_loops), std(logical_errors_loops), mean(failure_rate_loops), std(failure_rate_loops),
            mean(logical_errors_no_loops), std(logical_errors_no_loops), mean(failure_rate_no_loops), std(failure_rate_no_loops))
end

# L = 5 
# p = 0.1 
# max_loop_order = 0 
# numsamples = 10 
# pcmat = toric_code_X_parity_matrix(L)
# println(decode(pcmat, p, numsamples; pbias = 0.1, max_loop_order = 8))

# Parameters
Ls = [3,5,7]
ps = 0.00003:0.0003:0.003
numsamples = 5000
max_loop_order = 4

# Generate filename from command line arguments or default
filename_base = length(ARGS) > 0 ? ARGS[1] : "toric_code_results_$(Dates.format(now(), "yyyy-mm-dd_HH-MM-SS"))"
csv_filename = "data/$(filename_base).csv"


function red_gradient(n::Int)
    if n == 1
        return [RGB(1.0, 0.0, 0.0)]  # full red
    else
        return reverse([RGB(i, 0, 0) for i in range(0.4, stop=1.0, length=n)])
    end
end


# Initialize data collection
all_results = DataFrame(
    L = Int[],
    p = Float64[],
    logical_mean_loops = Float64[],
    logical_std_loops = Float64[],
    failure_mean_loops = Float64[],
    failure_std_loops = Float64[],
    logical_mean_no_loops = Float64[],
    logical_std_no_loops = Float64[],
    failure_mean_no_loops = Float64[],
    failure_std_no_loops = Float64[]
)

# Generate color gradient
colors = red_gradient(length(Ls))

# Create plots for both metrics
plt_logical = plot(title="Toric Code Logical Error Rate vs Physical Error Rate",
                  xlabel="Physical error rate p",
                  ylabel="Logical error rate",
                  legend=:topright)

plt_failure = plot(title="Toric Code Syndrome Failure Rate vs Physical Error Rate",
                  xlabel="Physical error rate p", 
                  ylabel="Syndrome failure rate",
                  legend=:topright)

# Loop over L values
for (i, L) in enumerate(Ls)
    println("Processing L = $L")
    pcmat = toric_code_X_parity_matrix(L)

    # With loop corrections
    logical_means_loops = Float64[]
    logical_stds_loops = Float64[]
    failure_means_loops = Float64[]
    failure_stds_loops = Float64[]
    
    # Without loop corrections
    logical_means_no_loops = Float64[]
    logical_stds_no_loops = Float64[]
    failure_means_no_loops = Float64[]
    failure_stds_no_loops = Float64[]

    # Create progress bar for this L value
    prog = Progress(length(ps), desc="L=$L: ")
    
    for (j, p) in enumerate(ps)
        # Decode both methods simultaneously
        μ_logical_loops, σ_logical_loops, μ_failure_loops, σ_failure_loops,
        μ_logical_no_loops, σ_logical_no_loops, μ_failure_no_loops, σ_failure_no_loops = decode(pcmat, p, numsamples, L; pbias=p, max_loop_order=max_loop_order)
        
        # Store results for loop method
        push!(logical_means_loops, μ_logical_loops)
        push!(logical_stds_loops, σ_logical_loops)
        push!(failure_means_loops, μ_failure_loops)
        push!(failure_stds_loops, σ_failure_loops)
        
        # Store results for no-loop method
        push!(logical_means_no_loops, μ_logical_no_loops)
        push!(logical_stds_no_loops, σ_logical_no_loops)
        push!(failure_means_no_loops, μ_failure_no_loops)
        push!(failure_stds_no_loops, σ_failure_no_loops)
        
        # Add to comprehensive data collection
        push!(all_results, (L, p, μ_logical_loops, σ_logical_loops, μ_failure_loops, σ_failure_loops,
                           μ_logical_no_loops, σ_logical_no_loops, μ_failure_no_loops, σ_failure_no_loops))
        
        # Update progress bar
        next!(prog, showvalues = [(:p, p), (:logical_loops, μ_logical_loops), (:logical_no_loops, μ_logical_no_loops),
                                  (:syndrome_loops, μ_failure_loops), (:syndrome_no_loops, μ_failure_no_loops)])
    end

    # Plot logical error rates
    plot!(plt_logical, ps, logical_means_loops;
          label="L = $L (with loops)",
          marker=:circle,
          color=colors[i],
          linestyle=:solid)
          
    plot!(plt_logical, ps, logical_means_no_loops;
          label="L = $L (no loops)",
          marker=:square,
          color=colors[i],
          linestyle=:dash)
          
    # Plot failure rates  
    plot!(plt_failure, ps, failure_means_loops;
          label="L = $L (with loops)",
          marker=:circle,
          color=colors[i],
          linestyle=:solid)
          
    plot!(plt_failure, ps, failure_means_no_loops;
          label="L = $L (no loops)",
          marker=:square,
          color=colors[i],
          linestyle=:dash)
end

# Save all results to CSV
CSV.write(csv_filename, all_results)
println("Results saved to $csv_filename")

# Create combined plot
plt_combined = plot(plt_logical, plt_failure, layout=(2,1), size=(800, 800))

# Save plots to visualization directory
mkpath("visualization")
savefig(plt_logical, "visualization/toric_code_logical_errors.png")
savefig(plt_failure, "visualization/toric_code_syndrome_failures.png") 
savefig(plt_combined, "visualization/toric_code_combined_metrics.png")

println("Plots saved to visualization/ directory")
println("Data saved to $csv_filename")
println("Press Enter to exit...")