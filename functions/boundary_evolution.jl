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
    return result
end

"""
    generate_ising_peps(Ny, Nx, beta, hz)

Generates a PEPS representing the partition function of a 2D classical Ising model.

The spins of the Ising model are located on the bonds of the PEPS lattice.
Each tensor at a vertex (y, x) is constructed by directly calculating the 
Boltzmann weight e^(-beta * H_local) for each spin configuration. H_local includes
diagonal nearest-neighbor ZZ interactions and a magnetic field `hz` term for each
connected spin.

Arguments:
- `Ny::Int`: The height of the lattice.
- `Nx::Int`: The width of the lattice.
- `beta::Float64`: The inverse temperature.
- `hz::Float64`: The external magnetic field strength.

Returns:
- A `Matrix{ITensor}` representing the Ising model PEPS.
"""
function generate_ising_peps(Ny::Int, Nx::Int, beta::Float64, hz::Float64)
    # Create the virtual indices, which represent the Ising spins on the bonds.
    # Index dimension 2 corresponds to spin up/down.
    h_inds = [Index(2, "Horz,y=$y,x=$x") for y in 1:Ny, x in 1:Nx-1]
    v_inds = [Index(2, "Vert,y=$y,x=$x") for y in 1:Ny-1, x in 1:Nx]

    peps = Matrix{ITensor}(undef, Ny, Nx)

    println("Generating Ising PEPS for a $Ny x $Nx lattice...")
    for y in 1:Ny
        for x in 1:Nx
            # --- Determine the indices for the current tensor ---
            u_ind = y > 1 ? v_inds[y-1, x] : nothing
            d_ind = y < Ny ? v_inds[y, x] : nothing
            l_ind = x > 1 ? h_inds[y, x-1] : nothing
            r_ind = x < Nx ? h_inds[y, x] : nothing

            valid_inds = filter(!isnothing, [u_ind, d_ind, l_ind, r_ind])
            
            if isempty(valid_inds) # Case for a 1x1 lattice
                peps[y,x] = ITensor(1.0) # No bonds, no spins, just a trivial tensor.
                continue
            end

            # Create the PEPS tensor for this site
            T = ITensor(valid_inds...)

            # Iterate over all possible spin configurations for this tensor.
            # CartesianIndices makes it easy to loop through all elements.
            for ci in CartesianIndices(T)
                # Helper to map an index value (1 or 2) to a spin value (+1 or -1).
                # Returns 0.0 if the index is not part of the current tensor (i.e., is `nothing`).
                function spin_val(ind)
                    if !isnothing(ind) && hasind(T, ind)
                        idxs = inds(T)
                        k = findfirst(x -> x == ind, idxs)
                        return ci[k] == 1 ? 1.0 : -1.0
                    else
                        return 0.0
                    end
                end

                su = spin_val(u_ind)
                sd = spin_val(d_ind)
                sl = spin_val(l_ind)
                sr = spin_val(r_ind)

                # --- Calculate local energy for this spin configuration ---
                # Field contribution: sum over all connected spins
                field_energy = -hz * (su + sd + sl + sr)

                # Interaction contribution: sum over diagonally adjacent pairs
                interaction_energy = -1.0 * (su*sl + su*sr + sd*sl + sd*sr)

                local_energy = field_energy + interaction_energy
                
                # Set the tensor element to the Boltzmann weight
                T[ci] = exp(-beta * local_energy)
            end
            peps[y,x] = T
        end
    end
    println("PEPS generation complete.")
    return peps
end

# """
#     generate_ising_peps(Ny, Nx, beta, hz)

# Generates a PEPS representing the partition function of a 2D classical Ising model.

# The spins of the Ising model are located on the bonds of the PEPS lattice.
# Each tensor at a vertex (y, x) encodes the Boltzmann weight e^(-beta * H_local),
# where H_local includes diagonal nearest-neighbor ZZ interactions and a magnetic
# field `hz` term for each connected spin.

# Arguments:
# - `Ny::Int`: The height of the lattice.
# - `Nx::Int`: The width of the lattice.
# - `beta::Float64`: The inverse temperature.
# - `hz::Float64`: The external magnetic field strength.

# Returns:
# - A `Matrix{ITensor}` representing the Ising model PEPS.
# """
# function generate_ising_peps(Ny::Int, Nx::Int, beta::Float64, hz::Float64)
#     # Create the virtual indices, which represent the Ising spins on the bonds.
#     h_inds = [Index(2, "Horz,y=$y,x=$x") for y in 1:Ny, x in 1:Nx-1]
#     v_inds = [Index(2, "Vert,y=$y,x=$x") for y in 1:Ny-1, x in 1:Nx]

#     peps = Matrix{ITensor}(undef, Ny, Nx)

#     println("Generating Ising PEPS for a $Ny x $Nx lattice...")
#     for y in 1:Ny
#         for x in 1:Nx
#             # --- Determine the indices and local Hamiltonian for each site ---
#             inds_list = Index[]
#             H_local = ITensor()

#             # Get the indices for the four potential directions (up, down, left, right)
#             u_ind = y > 1 ? v_inds[y-1, x] : nothing
#             d_ind = y < Ny ? v_inds[y, x] : nothing
#             l_ind = x > 1 ? h_inds[y, x-1] : nothing
#             r_ind = x < Nx ? h_inds[y, x] : nothing

#             # Collect the valid indices for the current tensor
#             local_inds_map = Dict("u" => u_ind, "d" => d_ind, "l" => l_ind, "r" => r_ind)
#             inds_list = [ind for ind in values(local_inds_map) if ind !== nothing]

#             if isempty(inds_list) # Case for a 1x1 lattice
#                 peps[y,x] = ITensor(1.0)
#                 continue
#             end
            
#             # --- Build the local Hamiltonian ---
#             # Initialize with the identity on the given indices
#             H_local = ITensor(inds_list...) 

#             # Add field terms: -hz * Z for each connected spin
#             for s_ind in inds_list
#                 # 
#                 H_local -= hz * pauliZ(s_ind)
#             end

#             # Add diagonal interaction terms: -1 * Z*Z
#             # This follows the model where the tensor encodes interactions
#             # between diagonally adjacent spins.
#             if !isnothing(u_ind) && !isnothing(l_ind); H_local -= pauliZ(u_ind) * pauliZ(l_ind); end
#             if !isnothing(u_ind) && !isnothing(r_ind); H_local -= pauliZ(u_ind) * pauliZ(r_ind); end
#             if !isnothing(d_ind) && !isnothing(l_ind); H_local -= pauliZ(d_ind) * pauliZ(l_ind); end
#             if !isnothing(d_ind) && !isnothing(r_ind); H_local -= pauliZ(d_ind) * pauliZ(r_ind); end

#             # The PEPS tensor is the Boltzmann weight e^(-beta * H)
#             # ITensor's exp function calculates the matrix exponential.
#             peps[y,x] = exp(-beta * H_local)
#         end
#     end
#     println("PEPS generation complete.")
#     return peps
# end

# # Define Pauli Z operator on a generic 2-level Index
# function pauliZ(ind::Index)
#     pind = prime(ind)
#     Z = ITensor(ind, pind)
#     Z[ind=>1, pind=>1] = 1
#     Z[ind=>2, pind=>2] = -1
#     return Z
# end


# # --- Main execution block ---
# function main()
#     # Define the dimensions of the PEPS
#     Ny = 3 # Height
#     Nx = 3 # Width
#     virt_dim = 2 # Virtual index dimension

#     # Create PEPS virtual indices
#     h_inds = [Index(virt_dim, "Horz,y=$y,x=$x") for y in 1:Ny, x in 1:Nx-1]
#     v_inds = [Index(virt_dim, "Vert,y=$y,x=$x") for y in 1:Ny-1, x in 1:Nx]

#     # Create a 2D array to hold the PEPS ITensors
#     peps = Matrix{ITensor}(undef, Ny, Nx)

#     # Fill the PEPS with random tensors (no physical indices)
#     println("Constructing random PEPS with no physical indices...")
#     for y in 1:Ny
#         for x in 1:Nx
#             inds_list = Index[]
#             if y > 1; push!(inds_list, dag(v_inds[y-1,x])); end # Top virtual link
#             if y < Ny; push!(inds_list, v_inds[y,x]); end # Bottom virtual link
#             if x > 1; push!(inds_list, dag(h_inds[y,x-1])); end # Left virtual link
#             if x < Nx; push!(inds_list, h_inds[y,x]); end # Right virtual link
            
#             if y == 1
#                 if x == 1
#                     peps[y,x] = ITensor(inds_list...)
#                     # TODO: generate tensor
#                     spin = i -> 3 - 2*i ## 1 -> 1, 2 -> -1
#                     for i1 in 1:2, i2 in 1:2, i3 in 1:2, i4 in 1:2
#                         s = spin(i1)*spin(i2) + spin(i2)*spin(i3) + spin(i3)*spin(i4) + spin(i4)*spin(i1) 
#                         s += h*(spin(i1) + spin(i2) + spin(i3) + spin(i4))/2 # divvy by two to avoid double counting
#                         peps[y,x][i1,i2,i3,i4] = exp(-β * s)
#                     end 
#                 end
#             end
#         end
#     end

#     # Perform the approximate contraction
#     result = contract_peps_no_phys(peps; cutoff=1E-6, maxdim=50)

#     println("\nApproximate PEPS Contraction Result: ", result)
# end
