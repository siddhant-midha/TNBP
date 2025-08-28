push!(LOAD_PATH, "../../functions/")
using BP
using Random, Plots, SparseArrays, ITensors, Statistics, ProgressMeter, Colors, LinearAlgebra, CSV, DataFrames, Dates, Graphs
include("../ldpc_tanner_loops.jl")
include("../../functions/tc_decode.jl")
include("../../functions/decoding.jl")


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


function toric_code_logical_operators(L::Int)
    # Compute the null space (logical operators) of the toric code parity check matrix
    pcmat = toric_code_X_parity_matrix(L)
    N = 2 * L^2
    logical_ops = zeros(Int, 2, N)
    
    # First logical operator: horizontal loop around the torus
    for i = 1:L
        logical_ops[1, (i-1)*L + 1] = 1  # First column of each row
    end
    
    # Second logical operator: vertical loop around the torus
    for j = 1:L
        logical_ops[2, L^2 + j] = 1  # First row of vertical edges
    end
    
    # Verify these are in the null space (PCM * logical_ops' = 0 mod 2)
    @assert all([(pcmat * logical_ops'[:,i]) .% 2 == zeros(L^2) for i in 1:2])
    
    return logical_ops
end



function interaction_tensor(disorder, β, indices)
    # Creates a Gibbs interaction tensor for spin variables
    # disorder: coupling value (J)
    # β: inverse temperature
    # indices: ITensor indices (each 2-dimensional representing spin ±1)
    
    it_tensor = ITensor(indices...)
    
    # Map from ITensor index values (1,2) to spin values (-1,1)
    map_value(x) = x == 1 ? -1 : x == 2 ? 1 : error("Unsupported index value: $x")
    
    # For each configuration in the tensor
    for I in eachindex(it_tensor)
        # Convert index values to spin values and compute energy
        spins = [map_value(I[i]) for i in 1:length(indices)]
        local_energy = disorder * prod(spins)  # J*s_i*s_j for 2 spins
        
        # Set Boltzmann weight
        it_tensor[I] = exp(β * local_energy)
    end 
    
    return it_tensor 
end


function periodic_square_lattice(L::Int)
    g = Graph(L^2)

    # Correct row-major indexing
    idx(i, j) = (i - 1) * L + j
    
    for i in 1:L
        for j in 1:L
            current = idx(i, j)
            right = j < L ? idx(i, j + 1) : idx(i, 1)  # Wrap rightmost to first column
            down = i < L ? idx(i + 1, j) : idx(1, j)    # Wrap bottom to top row
            add_edge!(g, current, right)
            add_edge!(g, current, down)
        end
    end
    return g
end

function get_partition_function_tensor(error_chain, L, β)
    # Construct toric code density matrix tensor network
    # L: linear size of the torus
    # error_chain: binary array representing errors (0 -> -1, 1 -> +1)
    # β: inverse temperature
    
    inv_idx(k) = ((k - 1) ÷ L + 1, (k - 1) % L + 1)
    idx(i, j) = (i - 1) * L + j
    
    # Function to map error (0,1) to disorder (-1,+1)
    error_to_disorder(e) = e == 0 ? -1.0 : 1.0

    g = periodic_square_lattice(L)
    N = L^2
    indmat = [Index(2, "i$(i)j$(j)") for i in 1:N, j in 1:N]
    site_tensors = []
    edge_tensors = []
    
    for index = 1:N 
        nbrs = neighbors(g, index)
        indices = [indmat[index, nbr] for nbr in nbrs]
        site_tensor = delta(indices)
        push!(site_tensors, site_tensor)
        
        ii, jj = inv_idx(index)
        
        # Handle horizontal edge (to the right)
        nbr1 = idx((ii) % L + 1, jj)
        # Get the corresponding error for horizontal edge (index is horizontal edge index)
        h_edge_idx = (ii - 1) * L + jj
        disorder = error_to_disorder(error_chain[h_edge_idx])
        push!(edge_tensors, interaction_tensor(disorder, β, [indmat[index, nbr1], indmat[nbr1, index]]))
        
        # Handle vertical edge (to the top)
        nbr2 = idx(ii, (jj) % L + 1)
        # Get the corresponding error for vertical edge (L^2 + index is vertical edge index)
        v_edge_idx = L^2 + (ii - 1) * L + jj
        disorder = error_to_disorder(error_chain[v_edge_idx])
        push!(edge_tensors, interaction_tensor(disorder, β, [indmat[index, nbr2], indmat[nbr2, index]]))
    end 

    return site_tensors, edge_tensors
end


function simple_match_syndrome(syndrome, L)
    # Simple function to find any error chain consistent with the syndrome
    # Connects pairs of excitations with straight paths (no optimization)
    
    # Initialize error configuration
    errors = zeros(Int, 2*L^2)
    
    # Find syndrome locations (vertex coordinates)
    excitations = []
    for i in 1:length(syndrome)
        if syndrome[i] == 1
            row = div(i-1, L) + 1
            col = mod1(i, L)
            push!(excitations, (row, col))
        end
    end
    
    # If no excitations, return zero error
    if isempty(excitations)
        return errors
    end
    
    # If odd number of excitations, this is impossible for a valid syndrome
    if length(excitations) % 2 != 0
        error("Invalid syndrome: odd number of excitations")
    end
    
    # Connect consecutive pairs of excitations with simple paths
    for i in 1:2:length(excitations)
        p1 = excitations[i]
        p2 = excitations[i+1]
        
        # First move horizontally from p1 to same column as p2
        r1, c1 = p1
        r2, c2 = p2
        
        # Horizontal path
        current_row = r1
        target_row = r2
        current_col = c1
        
        while current_row != target_row
            # Move one step toward target row
            if current_row < target_row
                next_row = current_row + 1
            else
                next_row = current_row - 1
            end
            
            # Handle wrapping
            if next_row == 0
                next_row = L
            elseif next_row == L + 1
                next_row = 1
            end
            
            # Add horizontal edge between current_row and next_row at current_col
            edge_row = min(current_row, next_row)
            if (current_row == L && next_row == 1) || (current_row == 1 && next_row == L)
                edge_row = L  # Handle wrap-around edge
            end
            
            h_edge_idx = (edge_row - 1) * L + current_col
            errors[h_edge_idx] = 1 - errors[h_edge_idx]
            
            current_row = next_row
        end
        
        # Then move vertically to reach p2
        while current_col != c2
            # Move one step toward target column
            if current_col < c2
                next_col = current_col + 1
            else
                next_col = current_col - 1
            end
            
            # Handle wrapping
            if next_col == 0
                next_col = L
            elseif next_col == L + 1
                next_col = 1
            end
            
            # Add vertical edge between current_col and next_col at current_row
            edge_col = min(current_col, next_col)
            if (current_col == L && next_col == 1) || (current_col == 1 && next_col == L)
                edge_col = L  # Handle wrap-around edge
            end
            
            v_edge_idx = L^2 + (current_row - 1) * L + edge_col
            errors[v_edge_idx] = 1 - errors[v_edge_idx]
            
            current_col = next_col
        end
    end
    
    return errors
end

function plot_error_chain(error_chain, L)
    """
    Plot an error chain on a toric code lattice.
    Blue edges indicate errors, black edges indicate no errors.
    """
    # Create figure
    p = plot(aspect_ratio=:equal, legend=false, grid=false, 
             xlims=(-0.5, L+0.5), ylims=(-0.5, L+0.5),
             title="Toric Code Error Chain (L=$L)")
    
    # Draw vertices
    for i in 1:L
        for j in 1:L
            scatter!([j], [i], color=:red, markersize=4)
        end
    end
    
    # Draw horizontal edges
    for i in 1:L
        for j in 1:L
            h_edge_idx = (i-1)*L + j
            color = error_chain[h_edge_idx] == 1 ? :blue : :black
            linewidth = error_chain[h_edge_idx] == 1 ? 3 : 1
            
            # Edge from (i,j) to (i+1,j) with periodic boundary
            next_i = i == L ? 1 : i + 1
            plot!([j, j], [i, next_i], color=color, linewidth=linewidth)
        end
    end
    
    # Draw vertical edges  
    for i in 1:L
        for j in 1:L
            v_edge_idx = L^2 + (i-1)*L + j
            color = error_chain[v_edge_idx] == 1 ? :blue : :black
            linewidth = error_chain[v_edge_idx] == 1 ? 3 : 1
            
            # Edge from (i,j) to (i,j+1) with periodic boundary
            next_j = j == L ? 1 : j + 1
            plot!([j, next_j], [i, i], color=color, linewidth=linewidth)
        end
    end
    
    return p
end


function plot_syndrome(syndrome, L)
    """
    Plot a syndrome on a toric code lattice.
    Red vertices indicate syndrome violations (excitations), black vertices indicate no violations.
    """
    # Create figure
    p = plot(aspect_ratio=:equal, legend=false, grid=false, 
             xlims=(-0.5, L+0.5), ylims=(-0.5, L+0.5),
             title="Toric Code Syndrome (L=$L)")
    
    # Draw all edges in light gray
    for i in 1:L
        for j in 1:L
            # Horizontal edges
            next_i = i == L ? 1 : i + 1
            plot!([j, j], [i, next_i], color=:lightgray, linewidth=1)
            
            # Vertical edges
            next_j = j == L ? 1 : j + 1
            plot!([j, next_j], [i, i], color=:lightgray, linewidth=1)
        end
    end
    
    # Draw vertices with syndrome information
    for i in 1:L
        for j in 1:L
            vertex_idx = (i-1)*L + j
            if syndrome[vertex_idx] == 1
                # Syndrome violation (excitation)
                scatter!([j], [i], color=:red, markersize=8, markerstrokewidth=2)
            else
                # No syndrome violation
                scatter!([j], [i], color=:black, markersize=4)
            end
        end
    end
    
    return p
end


function plot_syndrome_and_error(syndrome, error_chain, L)
    """
    Plot both syndrome and error chain on a toric code lattice.
    Red vertices indicate syndrome violations (excitations).
    Blue edges indicate errors, black edges indicate no errors.
    """
    # Create figure
    p = plot(aspect_ratio=:equal, legend=false, grid=false, 
             xlims=(-0.5, L+0.5), ylims=(-0.5, L+0.5),
             title="Toric Code: Syndrome + Error Chain (L=$L)")
    
    # Draw horizontal edges
    for i in 1:L
        for j in 1:L
            h_edge_idx = (i-1)*L + j
            color = error_chain[h_edge_idx] == 1 ? :blue : :lightgray
            linewidth = error_chain[h_edge_idx] == 1 ? 3 : 1
            
            # Edge from (i,j) to (i+1,j) with periodic boundary
            next_i = i == L ? 1 : i + 1
            plot!([j, j], [i, next_i], color=color, linewidth=linewidth)
        end
    end
    
    # Draw vertical edges  
    for i in 1:L
        for j in 1:L
            v_edge_idx = L^2 + (i-1)*L + j
            color = error_chain[v_edge_idx] == 1 ? :blue : :lightgray
            linewidth = error_chain[v_edge_idx] == 1 ? 3 : 1
            
            # Edge from (i,j) to (i,j+1) with periodic boundary
            next_j = j == L ? 1 : j + 1
            plot!([j, next_j], [i, i], color=color, linewidth=linewidth)
        end
    end
    
    # Draw vertices with syndrome information
    for i in 1:L
        for j in 1:L
            vertex_idx = (i-1)*L + j
            if syndrome[vertex_idx] == 1
                # Syndrome violation (excitation)
                scatter!([j], [i], color=:red, markersize=8, markerstrokewidth=2)
            else
                # No syndrome violation
                scatter!([j], [i], color=:black, markersize=4)
            end
        end
    end
    
    return p
end


function BP_free_energy(matched_errors, L, β)
    site_tensors, edge_tensors = get_partition_function_tensor(matched_errors, L, β)
    tensors = vcat(site_tensors, edge_tensors)
    adj_mat, edges, links = BP.get_adj_mat(tensors)
    messages = BP.get_messages(tensors, edges, links; random_part=0.01)
    messages = BP.message_passing(tensors, messages, edges, adj_mat; α=0.9, max_iters=500)
    Z_l = BP.get_fixed_point_list(tensors, messages, adj_mat)
    return sum(log.(real.(Z_l)))
end

function β_(p)
    ratio = p/(1.0 - p)
    return -1/2 * log(ratio)
end