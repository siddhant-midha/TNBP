using SparseArrays, Random, LinearAlgebra
include("../functions/tc_decode.jl")

# # ---------------- Your PCM (unchanged) ----------------
# function toric_code_X_parity_matrix(L::Int)
#     N = 2 * L^2   # edges: first L^2 horizontals, then L^2 verticals
#     M = L^2       # X checks (vertices)

#     pcmat = spzeros(Int, M, N)
#     for v in 1:M
#         i = div(v-1, L) + 1        # vertex row
#         j = mod1(v, L)             # vertex col
#         # horizontals touching (i,j): rows i-1 and i, same col j
#         h_edge1 = (mod1(i-1, L) - 1)*L + j
#         h_edge2 = (mod1(i,   L) - 1)*L + j
#         # verticals touching (i,j): row i, cols j-1 and j
#         v_edge1 = L^2 + (i-1)*L + mod1(j-1, L)
#         v_edge2 = L^2 + (i-1)*L + j
#         pcmat[v, h_edge1] = 1
#         pcmat[v, h_edge2] = 1
#         pcmat[v, v_edge1] = 1
#         pcmat[v, v_edge2] = 1
#     end
#     return pcmat
# end

# # group 1 (first L^2): H(row, col)   -- actually vertical edges
# @inline e_h(row::Int, col::Int, L::Int) = ((row-1) % L) * L + ((col-1) % L) + 1
# # group 2 (next  L^2): V(row, col)   -- actually horizontal edges
# @inline e_v(row::Int, col::Int, L::Int) = L^2 + ((row-1) % L) * L + ((col-1) % L) + 1

# # CORRECT plaquette boundary for PCM star at (i,j) (lower-left corner)
# plaquette_edges(i::Int, j::Int, L::Int) = [
#     e_h(i,            j,            L),  # left   vertical
#     e_h(i,            mod1(j+1, L), L),  # right  vertical
#     e_v(i,            j,            L),  # bottom horizontal
#     e_v(mod1(i+1, L), j,            L),  # top    horizontal
# ]

# # Build P: each row r=(i-1)*L+j has 1s on ∂P(i,j)
# function plaquette_boundary_matrix(L::Int)
#     rows = Int[]; cols = Int[]; vals = Int[]
#     for i in 1:L, j in 1:L
#         r = (i-1)*L + j
#         for e in plaquette_edges(i,j,L)
#             push!(rows, r); push!(cols, e); push!(vals, 1)
#         end
#     end
#     sparse(rows, cols, vals, L^2, 2L^2)
# end

# # Symmetric difference of two edge sets
# function symdiff(a::Vector{Int}, b::Vector{Int})
#     seen = Dict{Int,Int}()
#     @inbounds for v in a; seen[v] = get(seen,v,0) ⊻ 1; end
#     @inbounds for v in b; seen[v] = get(seen,v,0) ⊻ 1; end
#     [k for (k,bit) in seen if bit==1]
# end

# # Convert edge index set to 0/1 vector of length 2L^2
# edges_to_vec(S::Vector{Int}, L::Int) = begin
#     x = falses(2L^2); @inbounds for e in S; x[e] = .!x[e]; end; x
# end

# # --- GF(2) linear solve: does A * x == b have a solution? (Boolean) ---
# # Here we use a plain dense elimination on a copy of A (BitMatrix) for robustness.
# function gf2_solvable(A::SparseMatrixCSC{Int,Int}, b::BitVector)::Bool
#     m, n = size(A)
#     Ad = BitMatrix(A .% 2)         # m x n
#     bd = copy(b)                   # length m
#     r = 1                          # current row in elimination
#     for c in 1:n
#         # find pivot row at/after r with Ad[p,c]==1
#         p = 0
#         for i in r:m
#             if Ad[i,c]; p = i; break; end
#         end
#         if p == 0; continue; end   # no pivot in this column
#         # swap rows p and r
#         if p != r
#             Ad[p, :], Ad[r, :] = Ad[r, :], Ad[p, :]
#             bd[p], bd[r] = bd[r], bd[p]
#         end
#         # eliminate below
#         for i in r+1:m
#             if Ad[i,c]
#                 Ad[i, :] .⊻= Ad[r, :]
#                 bd[i] ⊻= bd[r]
#             end
#         end
#         r += 1
#         if r > m; break; end
#     end
#     # check consistency (0 == rhs?) for all-zero rows
#     for i in r:m
#         if all(!Ad[i,j] for j in 1:n) && bd[i]
#             return false
#         end
#     end
#     return true
# end

# # ---------------- Logical Z cycles as edge sets ----------------
# # H(j0): all horizontal edges at fixed column j0
# build_Hcycle_at_col(L::Int; j0::Int=1) = [e_h(i, ((j0-1)%L)+1, L) for i in 1:L]
# # V(i0): all vertical edges at fixed row i0
# build_Vcycle_at_row(L::Int; i0::Int=1) = [e_v(((i0-1)%L)+1, j, L) for j in 1:L]

# @inline function even_parity_on(S::Vector{Int}, support::Vector{Int})
#     supp = BitSet(support); p = 0
#     @inbounds for e in S; if e in supp; p ⊻= 1; end; end
#     p == 0
# end

# # ---------------- Checks ----------------
# """
#     check_decode_edges(errors, errors_true, L; j0=1, i0=1)

# Returns (syndrome_zero::Bool, homology_trivial::Bool) for Sdiff = errors ⊻ errors_true,
# using your PCM for syndrome and cycle parities for homology.
# """
# # function check_decode_edges(errors::Vector{Int}, errors_true::Vector{Int}, L::Int; j0::Int=1, i0::Int=1)
# #     pcmat = toric_code_X_parity_matrix(L)
# #     Sdiff = symdiff(errors, errors_true)

# #     s = edges_to_vec(Sdiff, L)             # 0/1 vector length 2L^2
# #     syn = pcmat * Int.(s)                  # integer vector; reduce mod 2
# #     syndrome_zero = all(x -> (x & 1) == 0, syn)

# #     H = edges_to_vec(build_Hcycle_at_col(L; j0=j0), L)
# #     V = edges_to_vec(build_Vcycle_at_row(L; i0=i0), L)
# #     wy_zero = (sum(Int.(H) .* Int.(s)) & 1) == 0   # even parity vs H(j0)
# #     wx_zero = (sum(Int.(V) .* Int.(s)) & 1) == 0   # even parity vs V(i0)
# #     homology_trivial = wx_zero && wy_zero

# #     return syndrome_zero, homology_trivial
# # end
# # --- main check using PCM for syndrome and P for homology ---
# function check_decode_edges(errors::Vector{Int}, errors_true::Vector{Int}, L::Int)
#     pcmat = toric_code_X_parity_matrix(L)
#     P = plaquette_boundary_matrix(L)

#     Sdiff = symdiff(errors, errors_true)
#     s = edges_to_vec(Sdiff, L)
#     syn = pcmat * Int.(s)
#     syndrome_zero = all(x -> (x & 1) == 0, syn)

#     # homology trivial iff Sdiff is a sum of plaquette boundaries: Pᵀ x = s has a solution
#     homology_trivial = gf2_solvable(sparse(transpose(P)), BitVector(s))  # A = Pᵀ (size 2L^2 × L^2)

#     return syndrome_zero, homology_trivial
# end

# # ---------------- Plaquette boundary (for random test) ----------------
# function plaquette_edges(i::Int, j::Int, L::Int)
#     i0 = ((i-1)%L)+1; j0 = ((j-1)%L)+1
#     h_bot = e_h(((i0-1-1)%L)+1, j0, L)  # row i-1, col j
#     h_top = e_h(i0, j0, L)              # row i,   col j
#     v_left  = e_v(i0, ((j0-1-1)%L)+1, L)# row i,   col j-1
#     v_right = e_v(i0, j0, L)            # row i,   col j
#     [h_bot, h_top, v_left, v_right]
# end

bad_vertices(S::Vector{Int}, L::Int) = begin
    pcmat = toric_code_X_parity_matrix(L)
    s = edges_to_vec(S, L)
    syn = pcmat * Int.(s)
    [(div(v-1,L)+1, mod1(v,L)) for v in findall(x->(x & 1)==1, syn)]
end

# random edge set
function random_errors_true_edges(L::Int; p=0.15, rng=Random.default_rng())
    S = Int[]
    @inbounds for i in 1:L, j in 1:L
        if rand(rng) < p; push!(S, e_h(i,j,L)); end
        if rand(rng) < p; push!(S, e_v(i,j,L)); end
    end
    S
end

# remove one edge from a loop to make it open
make_incomplete(edges::Vector{Int}) = edges[1:end-1]

# ---- RUN TESTS ----
function run_random_homology_tests(L::Int=9; rng=Random.MersenneTwister(123))
    errors_true = random_errors_true_edges(L; p=0.20, rng)

    # pick reference cycles
    j0 = rand(rng, 1:L)
    i0 = rand(rng, 1:L)
    H0 = build_Hcycle_at_col(L; j0=j0)     # non-contractible (y-winding)
    V0 = build_Vcycle_at_row(L; i0=i0)     # non-contractible (x-winding)

    println("L=$L, j0=$j0 (H-cycle@col), i0=$i0 (V-cycle@row)")
    println("sizes: |errors_true|=$(length(errors_true)), |H0|=$(length(H0)), |V0|=$(length(V0))")

    # (1) add two loops (same loop twice) -> cancels mod 2 → (true, true)
    e1 = symdiff(errors_true, H0)                    # add H0
    e1 = symdiff(e1, H0)                             # add H0 again (cancels)
    println("(1) add same loop twice: ", check_decode_edges(e1, errors_true, L))

    # (1b) add two parallel disjoint loops H(j0) ⊕ H(j1) -> trivial homology → (true, true)
    j1 = (j0 % L) + 1
    H1 = build_Hcycle_at_col(L; j0=j1)
    e1b = symdiff(errors_true, symdiff(H0, H1))
    println("(1b) add two parallel loops H(j0),H(j1): ", check_decode_edges(e1b, errors_true, L), " expected (true, true)")

    # (2) add two disjoint loops of different orientation H ⊕ V -> change class → (true, false)
    e2 = symdiff(errors_true, symdiff(H0, V0))
    println("(2) add H ⊕ V: ", check_decode_edges(e2, errors_true, L), "expected (true, false)")

    # (3) add one incomplete loop (open chain) -> nonzero syndrome → (false, false)
    incH = make_incomplete(H0)   # remove one edge from H0
    e3 = symdiff(errors_true, incH)
    println("(3) add incomplete H: ", check_decode_edges(e3, errors_true, L), "expected (false, false)")

    # Optionally verify that a single plaquette still yields (true,true)
    iP, jP = rand(rng, 1:L), rand(rng, 1:L)
    loopP = plaquette_edges(iP, jP, L)
    eP = symdiff(errors_true, loopP)
    println("(P) add 1 plaquette: ", check_decode_edges(eP, errors_true, L), "expected (true, true)")
end

# 2-plaquette loop (rectangle of size 1×2 plaquettes)
function two_plaquette_loop(i::Int, j::Int, L::Int; orientation::Symbol=:horizontal)
    if orientation === :horizontal
        # plaquettes at (i,j) and (i, j+1)
        return symdiff(plaquette_edges(i, j, L), plaquette_edges(i, mod1(j+1, L), L))
    elseif orientation === :vertical
        # plaquettes at (i,j) and (i+1, j)
        return symdiff(plaquette_edges(i, j, L), plaquette_edges(mod1(i+1, L), j, L))
    else
        error("orientation must be :horizontal or :vertical")
    end
end

# Optional: quick checker to confirm it's a 6-edge loop with zero-star syndrome
bad_vertices(S::Vector{Int}, L::Int) = begin
    pcmat = toric_code_X_parity_matrix(L)
    x = falses(2L^2); for e in S; x[e] = .!x[e]; end
    syn = pcmat * Int.(x)
    [(div(v-1,L)+1, mod1(v,L)) for v in findall(y->(y & 1)==1, syn)]
end

# ---- Test: add a big loop made of two plaquettes ----
function test_two_plaquette_loop(L::Int=9, seed::Int=123)
    rng = Random.MersenneTwister(seed)
    # random baseline errors_true (edges)
    function random_edges(L; p=0.2, rng=Random.default_rng())
        S = Int[]
        for i in 1:L, j in 1:L
            if rand(rng) < p; push!(S, e_h(i,j,L)); end
            if rand(rng) < p; push!(S, e_v(i,j,L)); end
        end
        S
    end
    errors_true = random_edges(L; p=0.20, rng)

    i, j = rand(rng, 1:L), rand(rng, 1:L)

    # horizontal 1×2 loop
    loopH = two_plaquette_loop(i, j, L; orientation=:horizontal)
    println("Horizontal 1×2 loop:")
    println("  length(loopH) = $(length(loopH)), expected = 6")
    println("  bad_vertices(loopH, L) = $(bad_vertices(loopH, L)), expected = []")
    eH = symdiff(errors_true, loopH)
    resultH = check_decode_edges(eH, errors_true, L)
    println("  check_decode_edges(eH, errors_true, L) = $resultH, expected = (true, true)")

    # vertical 2×1 loop
    loopV = two_plaquette_loop(i, j, L; orientation=:vertical)
    println("Vertical 2×1 loop:")
    println("  length(loopV) = $(length(loopV)), expected = 6")
    println("  bad_vertices(loopV, L) = $(bad_vertices(loopV, L)), expected = []")
    eV = symdiff(errors_true, loopV)
    resultV = check_decode_edges(eV, errors_true, L)
    println("  check_decode_edges(eV, errors_true, L) = $resultV, expected = (true, true)")
end


# ---------------- Tests ----------------
if abspath(PROGRAM_FILE) == @__FILE__
    L = 5
    errors_true = Int[]
    errors_dec  = Int[]

    # (true, true) — empty
    @show check_decode_edges(errors_dec, errors_true, L)

    # (true, false) — a closed nontrivial cycle:
    # choose H(j0): all horizontals in one column (commutes with all stars; flips H parity)
    loop_nontrivial = build_Hcycle_at_col(L; j0=1)
    @show check_decode_edges(loop_nontrivial, errors_true, L)

    # (false, true) — open local pattern that hits neither cycle oddly
    open_trivial = [e_v(2,3,L), e_v(3,3,L)]
    @show check_decode_edges(open_trivial, errors_true, L)

    # (false, false) — open pattern that flips at least one cycle parity:
    # take a single edge from H(j0) to make H-parity odd, plus one unrelated edge to keep it open
    open_nontrivial = [e_h(1,1,L), e_v(2,3,L)]
    @show check_decode_edges(open_nontrivial, errors_true, L)

    # Random plaquette boundary keeps (true, true)
    L = 7
    rng = MersenneTwister(42)

    # random edge set
    function random_edges(L; p=0.2, rng=Random.default_rng())
        S = Int[]
        for i in 1:L, j in 1:L
            if rand(rng) < p; push!(S, e_h(i,j,L)); end
            if rand(rng) < p; push!(S, e_v(i,j,L)); end
        end
        S
    end

    errors_true = random_edges(L; p=0.20, rng)
    i = rand(rng, 1:L); j = rand(rng, 1:L)
    loop = plaquette_edges(i, j, L)
    errors_dec = symdiff(errors_true, loop)
    Sdiff = symdiff(errors_dec, errors_true)

    println("Sdiff == plaquette? ", Set(Sdiff) == Set(loop))
    println("bad vertices: ", bad_vertices(Sdiff, L))
    println("check_decode_edges(errors_dec, errors_true, L) => ",
        check_decode_edges(errors_dec, errors_true, L))   # -> (true, true)

    run_random_homology_tests()
    test_two_plaquette_loop()
end
