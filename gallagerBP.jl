using Serialization
using ITensors, ITensorMPS, Plots
using ProgressMeter, Statistics, Random, SparseArrays
using JLD2, Printf
include("ldpc_tanner_loops.jl")
push!(LOAD_PATH, "functions/")
using BP

function tensorargmax(probs)
    ix = inds(probs)[1]
    probs = [real((probs)[ix=>n]) for n in 1:dim(ix)]
    probs ./= sum(probs)
    return Int(argmax(probs)-1) #sample_bit(probs[1])
end 


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

"""
    gallager_ldpc_matrix(n::Int, d_v::Int, d_c::Int)

Generate a regular (d_v, d_c) LDPC parity-check matrix using Gallager construction.

# Arguments
- `n`: Number of bits (columns)
- `d_v`: Number of ones per column (column weight)
- `d_c`: Number of ones per row (row weight)

# Returns
- `H`: A sparse parity-check matrix of size m × n
"""
function gallager_ldpc_matrix(n::Int, d_v::Int, d_c::Int)
    # Ensure total ones are compatible
    if n * d_v % d_c != 0
        error("n*d_v must be divisible by d_c")
    end
    m = (n * d_v) ÷ d_c  # number of rows

    # Initialize empty matrix
    H = spzeros(Int8, m, n)

    # First, build the base structure
    rows_per_group = m ÷ d_v
    for i in 0:(d_v-1)
        perm = randperm(n)
        for j in 1:rows_per_group
            row = i * rows_per_group + j
            cols = perm[((j-1)*d_c + 1):(j*d_c)]
            H[row, cols] .= 1
        end
    end

    return H
end
 
function get_marginal(tensors,adj_mat,messages,index)
    nbrs = BP.get_nbrs(adj_mat, index)
    Z_local = tensors[index] 
    for nbr in nbrs
        Z_local *= messages[nbr,index] 
    end
    return Z_local
end 


function loops_regular_code_sim_batched(n, d_v, d_c, p, num_batches, samples_per_batch; pbias = 0.1, max_loop_order = 4)
    ## keeps a list of priors on all bits, and adds to the prior as and when a loop is encountered
    ## should work for arbitrary loop orders
    println("p = $p, N = $n, batches = $num_batches, samples per batch = $samples_per_batch")
    batch_error_rates_loops = zeros(num_batches)
    batch_error_rates_no_loops = zeros(num_batches)
    for batch in 1:num_batches
        println("batch=$batch")
        pcmat = gallager_ldpc_matrix(n, d_v, d_c)
        tannerloopslist = [find_tanner_loops(pcmat, d; max_length=max_loop_order) for d in 1:n]
        # Track errors within this batch
        batch_logical_errors_loops = 0
        batch_logical_errors_no_loops = 0
        
        for samp in 1:samples_per_batch
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
            messages = BP.message_passing(tensors,messages,edges,adj_mat;α=0.95, max_iters=500,diagnose=false,normalise=true)
            
            # Decode each data qubit
            for d = 1:n 
                probs = get_marginal(vcat(data_tensors, syn_tensors), adj_mat, messages, d)
                vacuum = tensorargmax(probs)
                errors_no_loops[d] = vacuum

                # Get loops for this data qubit
                tannerloops = tannerloopslist[d]
                loop_list = [tannerloop.edges for tannerloop in tannerloops] 
                data_bits_involved_list = [tannerloop.data_bits for tannerloop in tannerloops]
                check_bits_involved_list = [tannerloop.check_bits for tannerloop in tannerloops]
                
                loopprobs = ITensor([0,0], data_indices[d])
                for (i, loop) in enumerate(loop_list)
                    data_bits_involved = data_bits_involved_list[i] 
                    check_bits_involved = check_bits_involved_list[i] 
                    mtensors = vcat(get_marginal_data_tensors(data_tensors, data_indices, errors_loops; exclude=[d]), syn_tensors)
                    
                    if !isempty(setdiff(data_bits_involved, [d]))
                        normlz1 = scalar(prod([get_marginal(mtensors, adj_mat, messages, other_data_bit) 
                                            for other_data_bit in collect(setdiff(data_bits_involved, [d]))]))
                    else
                        normlz1 = 1.0
                    end
                    
                    if !isempty(check_bits_involved)
                        normlz2 = scalar(prod([get_marginal(mtensors, adj_mat, messages, n+check_bit) 
                                            for check_bit in check_bits_involved]))
                    else
                        normlz2 = 1.0
                    end
                    
                    normlz = normlz1 * normlz2
                    if normlz != 0
                        change = BP.loop_contribution(loop, messages, mtensors, edges, links, adj_mat) / normlz
                        loopprobs += change 
                    end
                end
                
                loopcorr = tensorargmax(probs + loopprobs)
                errors_loops[d] = loopcorr

                # Early stopping: if both decoders made errors, no point continuing
                if (errors_no_loops[d] != errors_true[d]) && (errors_loops[d] != errors_true[d]) 
                    break
                end 

                
            end
            
            batch_logical_errors_loops += (sum(errors_loops .!= errors_true) > 0)
            batch_logical_errors_no_loops += (sum(errors_no_loops .!= errors_true) > 0)
        end
        
        # Store error rates for this batch
        batch_error_rates_loops[batch] = batch_logical_errors_loops / samples_per_batch
        batch_error_rates_no_loops[batch] = batch_logical_errors_no_loops / samples_per_batch
    end
    
    # Return overall average error rates across all batches
    return (mean(batch_error_rates_loops), std(batch_error_rates_loops)), 
           (mean(batch_error_rates_no_loops), std(batch_error_rates_no_loops))
end



using Dates

# Get current timestamp for labeling
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

# Parameters
d_v = 3 
d_c = 4
num_batches = 50
samples_per_batch = 50
max_loop_order = 6
n_ar = [48, 96]
ps = LinRange(0.1, 0.2, 20)

# Get task parameters from SLURM array
task_id = parse(Int, ARGS[1])
p = ps[task_id]

@info "Task $task_id: Computing for p = $(p)"

# Run simulations for all n values
for n in n_ar
    @info "Computing for n = $(n), p = $(p)"
    
    # Run simulation
    (mean_loops, std_loops), (mean_no_loops, std_no_loops) = loops_regular_code_sim_batched(
        n, d_v, d_c, p, num_batches, samples_per_batch; 
        pbias = p, max_loop_order = max_loop_order
    )
    
    # Save results with detailed labeling
    filename = "gallagerBP0results/results_task$(task_id)_n$(n)_p$(round(p, digits=4))_$(timestamp).jld2"
    
    save(filename, 
         # Simulation parameters
         "n", n,
         "p", p,
         "d_v", d_v,
         "d_c", d_c,
         "num_batches", num_batches,
         "samples_per_batch", samples_per_batch,
         "max_loop_order", max_loop_order,
         "pbias", p,
         "task_id", task_id,
         "timestamp", timestamp,
         
         # Results with loops
         "mean_error_rate_loops", mean_loops,
         "std_error_rate_loops", std_loops,
         
         # Results without loops (vacuum)
         "mean_error_rate_no_loops", mean_no_loops,
         "std_error_rate_no_loops", std_no_loops
    )
    
    @info "Saved results to $filename"
    @info "  Loops decoder: $(mean_loops) ± $(std_loops)"
    @info "  No-loops decoder: $(mean_no_loops) ± $(std_no_loops)"
end

@info "Task $task_id completed successfully"

