push!(LOAD_PATH, "../functions/")
using BP
using Random, Plots, SparseArrays, ITensors, Statistics, ProgressMeter, Colors, LinearAlgebra

include("../tc.jl")
include("../ldpc_tanner_loops.jl")
include("../functions/tc_decode.jl")

function create_plaquette_error(L::Int, plaquette_i::Int, plaquette_j::Int)
    """
    Create an error pattern by flipping the left and top edges of a specific plaquette.
    For L=3 toric code with periodic boundaries.
    
    Plaquettes are indexed by (i,j) where i,j ∈ {1,2,3}
    Each plaquette has 4 edges: top, right, bottom, left
    We flip the left and top edges.
    """
    N = 2 * L^2  # Total number of qubits (edges)
    errors = zeros(Int, N)
    
    # For plaquette at (i,j), find the left and top edge indices
    # Top edge: horizontal edge above the plaquette 
    top_edge = (plaquette_i - 1) * L + plaquette_j
    
    # Left edge: vertical edge to the left of the plaquette
    left_edge = L^2 + (plaquette_i - 1) * L + mod1(plaquette_j - 1, L)
    
    # Flip these edges
    errors[top_edge] = 1
    errors[left_edge] = 1
    
    return errors
end

function print_priors(priors, data_indices, label)
    println("\n=== $label ===")
    for (i, prior) in enumerate(priors)
        ix = inds(prior)[1]
        probs = [real(prior[ix=>n]) for n in 1:dim(ix)]
        probs ./= sum(probs)
        println("Qubit $i: P(0)=$(round(probs[1], digits=6)), P(1)=$(round(probs[2], digits=6))")
    end
end

function test_toric_degeneracy()
    println("=== Toric Code Degeneracy Test ===")
    println("L=3 toric code, flipping left and top edges of plaquette (2,2)")
    
    L = 3
    pcmat = toric_code_X_parity_matrix(L)
    
    # Create error pattern: flip left and top edges of plaquette (2,2)
    errors_true = create_plaquette_error(L, 2, 2)
    
    println("\nError pattern (first 9 are horizontal edges, next 9 are vertical edges):")
    println("Horizontal edges: ", errors_true[1:L^2])
    println("Vertical edges: ", errors_true[L^2+1:2*L^2])
    
    # Compute syndrome
    syndrome = pcmat * errors_true .% 2
    println("\nSyndrome: ", syndrome)
    
    # Set up decoding problem
    pbias = 0.1
    data_tensors, syn_tensors, data_indices = get_network(pcmat, syndrome, pbias)
    
    m, n = size(pcmat)
    errors_loops = Int.(-1 .* ones(n))
    errors_no_loops = Int.(-1 .* ones(n))
    
    # Set up tensor network and BP
    tensors = vcat(get_marginal_data_tensors(data_tensors, data_indices, errors_loops), syn_tensors)
    adj_mat, edges, links = BP.get_adj_mat(tensors)
    messages = BP.get_messages(tensors, edges, links)
    messages = BP.message_passing(tensors, messages, edges, links, adj_mat; α=0.95, max_iters=500, diagnose=false, normalise=true)
    
    # Initialize priors for both methods
    priors_loops = [ITensor([0.,0.], data_indices[i]) for i in 1:n]
    priors_no_loops = [ITensor([0.,0.], data_indices[i]) for i in 1:n]
    
    # Get Tanner loops for loop corrections
    tannerloopslist = [find_tanner_loops(pcmat, d; max_length=8) for d in 1:n]
    
    println("\n=== Computing Priors ===")
    
    # Compute vacuum contributions (same for both methods)
    for d = 1:n
        probs = get_marginal(vcat(data_tensors, syn_tensors), adj_mat, messages, d)
        priors_loops[d] += probs
        priors_no_loops[d] += probs
    end
    
    print_priors(priors_no_loops, data_indices, "BP Only (Vacuum Contribution)")
    
    # Apply loop corrections
    for d = 1:n
        tannerloops = tannerloopslist[d]
        loop_list = [tannerloop.edges for tannerloop in tannerloops]
        data_bits_involved_list = [tannerloop.data_bits for tannerloop in tannerloops]
        check_bits_involved_list = [tannerloop.check_bits for tannerloop in tannerloops]
        
        for (i, loop) in enumerate(loop_list)
            data_bits_involved = data_bits_involved_list[i]
            check_bits_involved = check_bits_involved_list[i]
            data_bits_involved_other = collect(setdiff(data_bits_involved, [d]))
            
            for data_bit in data_bits_involved_other
                mtensors = vcat(get_marginal_data_tensors(data_tensors, data_indices, errors_loops; exclude=[data_bit]), syn_tensors)
                
                # Normalization from other nodes in the loop
                normlz = prod([get_marginal(mtensors, adj_mat, messages, other_data_bit) for other_data_bit in collect(setdiff(data_bits_involved, [data_bit]))])
                normlz *= prod([get_marginal(mtensors, adj_mat, messages, n + check_bit) for check_bit in check_bits_involved])
                
                # Update prior with loop contribution
                priors_loops[data_bit] += loop_contribution(loop, messages, mtensors, edges, links, adj_mat) / normlz
            end
        end
    end
    
    print_priors(priors_loops, data_indices, "BP + Loop Corrections")
    
    # Make hard decisions and compare
    println("\n=== Hard Decisions ===")
    
    decoded_bp = Int[]
    decoded_loops = Int[]
    
    for d = 1:n
        # BP only decision
        ix_bp = inds(priors_no_loops[d])[1]
        probs_bp = [real(priors_no_loops[d][ix_bp=>k]) for k in 1:dim(ix_bp)]
        probs_bp ./= sum(probs_bp)
        decision_bp = probs_bp[1] > 0.5 ? 0 : 1
        push!(decoded_bp, decision_bp)
        
        # BP + loops decision
        ix_loops = inds(priors_loops[d])[1]
        probs_loops = [real(priors_loops[d][ix_loops=>k]) for k in 1:dim(ix_loops)]
        probs_loops ./= sum(probs_loops)
        decision_loops = probs_loops[1] > 0.5 ? 0 : 1
        push!(decoded_loops, decision_loops)
    end
    
    println("True errors:     ", errors_true)
    println("BP decoded:      ", decoded_bp)
    println("BP+loops decoded:", decoded_loops)
    
    # Check syndrome satisfaction
    syndrome_bp = pcmat * decoded_bp .% 2
    syndrome_loops = pcmat * decoded_loops .% 2
    
    println("\nSyndrome check:")
    println("Target syndrome: ", syndrome)
    println("BP syndrome:     ", syndrome_bp)
    println("Loops syndrome:  ", syndrome_loops)
    
    println("BP syndrome match: ", all(syndrome_bp .== syndrome))
    println("Loops syndrome match: ", all(syndrome_loops .== syndrome))
    
    # Check logical error (homology)
    bp_success, bp_logical = check_decode_edges(decoded_bp, errors_true, L)
    loops_success, loops_logical = check_decode_edges(decoded_loops, errors_true, L)
    
    println("\nLogical error check:")
    println("BP logical success: ", bp_logical)
    println("Loops logical success: ", loops_logical)
end

# Run the test
test_toric_degeneracy()