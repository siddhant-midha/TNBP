using SparseArrays, Random, LinearAlgebra

# ---------------- Your PCM (unchanged) ----------------
function toric_code_X_parity_matrix(L::Int)
    N = 2 * L^2   # edges: first L^2 horizontals, then L^2 verticals
    M = L^2       # X checks (vertices)

    pcmat = spzeros(Int, M, N)
    for v in 1:M
        i = div(v-1, L) + 1        # vertex row
        j = mod1(v, L)             # vertex col
        # horizontals touching (i,j): rows i-1 and i, same col j
        h_edge1 = (mod1(i-1, L) - 1)*L + j
        h_edge2 = (mod1(i,   L) - 1)*L + j
        # verticals touching (i,j): row i, cols j-1 and j
        v_edge1 = L^2 + (i-1)*L + mod1(j-1, L)
        v_edge2 = L^2 + (i-1)*L + j
        pcmat[v, h_edge1] = 1
        pcmat[v, h_edge2] = 1
        pcmat[v, v_edge1] = 1
        pcmat[v, v_edge2] = 1
    end
    return pcmat
end

# group 1 (first L^2): H(row, col)   -- actually vertical edges
@inline e_h(row::Int, col::Int, L::Int) = ((row-1) % L) * L + ((col-1) % L) + 1
# group 2 (next  L^2): V(row, col)   -- actually horizontal edges
@inline e_v(row::Int, col::Int, L::Int) = L^2 + ((row-1) % L) * L + ((col-1) % L) + 1

# CORRECT plaquette boundary for PCM star at (i,j) (lower-left corner)
plaquette_edges(i::Int, j::Int, L::Int) = [
    e_h(i,            j,            L),  # left   vertical
    e_h(i,            mod1(j+1, L), L),  # right  vertical
    e_v(i,            j,            L),  # bottom horizontal
    e_v(mod1(i+1, L), j,            L),  # top    horizontal
]

# Build P: each row r=(i-1)*L+j has 1s on ∂P(i,j)
function plaquette_boundary_matrix(L::Int)
    rows = Int[]; cols = Int[]; vals = Int[]
    for i in 1:L, j in 1:L
        r = (i-1)*L + j
        for e in plaquette_edges(i,j,L)
            push!(rows, r); push!(cols, e); push!(vals, 1)
        end
    end
    sparse(rows, cols, vals, L^2, 2L^2)
end

# Symmetric difference of two edge sets
function symdiff(a::Vector{Int}, b::Vector{Int})
    seen = Dict{Int,Int}()
    @inbounds for v in a; seen[v] = get(seen,v,0) ⊻ 1; end
    @inbounds for v in b; seen[v] = get(seen,v,0) ⊻ 1; end
    [k for (k,bit) in seen if bit==1]
end

# Convert edge index set to 0/1 vector of length 2L^2
edges_to_vec(S::Vector{Int}, L::Int) = begin
    x = falses(2L^2); @inbounds for e in S; x[e] = .!x[e]; end; x
end

# --- GF(2) linear solve: does A * x == b have a solution? (Boolean) ---
# Here we use a plain dense elimination on a copy of A (BitMatrix) for robustness.
function gf2_solvable(A::SparseMatrixCSC{Int,Int}, b::BitVector)::Bool
    m, n = size(A)
    Ad = BitMatrix(A .% 2)         # m x n
    bd = copy(b)                   # length m
    r = 1                          # current row in elimination
    for c in 1:n
        # find pivot row at/after r with Ad[p,c]==1
        p = 0
        for i in r:m
            if Ad[i,c]; p = i; break; end
        end
        if p == 0; continue; end   # no pivot in this column
        # swap rows p and r
        if p != r
            Ad[p, :], Ad[r, :] = Ad[r, :], Ad[p, :]
            bd[p], bd[r] = bd[r], bd[p]
        end
        # eliminate below
        for i in r+1:m
            if Ad[i,c]
                Ad[i, :] .⊻= Ad[r, :]
                bd[i] ⊻= bd[r]
            end
        end
        r += 1
        if r > m; break; end
    end
    # check consistency (0 == rhs?) for all-zero rows
    for i in r:m
        if all(!Ad[i,j] for j in 1:n) && bd[i]
            return false
        end
    end
    return true
end

# ---------------- Logical Z cycles as edge sets ----------------
# H(j0): all horizontal edges at fixed column j0
build_Hcycle_at_col(L::Int; j0::Int=1) = [e_h(i, ((j0-1)%L)+1, L) for i in 1:L]
# V(i0): all vertical edges at fixed row i0
build_Vcycle_at_row(L::Int; i0::Int=1) = [e_v(((i0-1)%L)+1, j, L) for j in 1:L]

@inline function even_parity_on(S::Vector{Int}, support::Vector{Int})
    supp = BitSet(support); p = 0
    @inbounds for e in S; if e in supp; p ⊻= 1; end; end
    p == 0
end

# ---------------- Checks ----------------
"""
    check_decode_edges(errors, errors_true, L; j0=1, i0=1)

Returns (syndrome_zero::Bool, homology_trivial::Bool) for Sdiff = errors ⊻ errors_true,
using your PCM for syndrome and cycle parities for homology.
"""
# function check_decode_edges(errors::Vector{Int}, errors_true::Vector{Int}, L::Int; j0::Int=1, i0::Int=1)
#     pcmat = toric_code_X_parity_matrix(L)
#     Sdiff = symdiff(errors, errors_true)

#     s = edges_to_vec(Sdiff, L)             # 0/1 vector length 2L^2
#     syn = pcmat * Int.(s)                  # integer vector; reduce mod 2
#     syndrome_zero = all(x -> (x & 1) == 0, syn)

#     H = edges_to_vec(build_Hcycle_at_col(L; j0=j0), L)
#     V = edges_to_vec(build_Vcycle_at_row(L; i0=i0), L)
#     wy_zero = (sum(Int.(H) .* Int.(s)) & 1) == 0   # even parity vs H(j0)
#     wx_zero = (sum(Int.(V) .* Int.(s)) & 1) == 0   # even parity vs V(i0)
#     homology_trivial = wx_zero && wy_zero

#     return syndrome_zero, homology_trivial
# end
# --- main check using PCM for syndrome and P for homology ---
function check_decode_edges(errors::Vector{Int}, errors_true::Vector{Int}, L::Int)
    pcmat = toric_code_X_parity_matrix(L)
    P = plaquette_boundary_matrix(L)

    Sdiff = symdiff(errors, errors_true)
    s = edges_to_vec(Sdiff, L)
    syn = pcmat * Int.(s)
    syndrome_zero = all(x -> (x & 1) == 0, syn)

    # homology trivial iff Sdiff is a sum of plaquette boundaries: Pᵀ x = s has a solution
    homology_trivial = gf2_solvable(sparse(transpose(P)), BitVector(s))  # A = Pᵀ (size 2L^2 × L^2)

    return syndrome_zero, homology_trivial
end