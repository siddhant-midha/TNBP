# ---------- indexing helpers ----------
@inline idx(x::Int, y::Int, L::Int) = (x-1) % L * L + (y-1) % L + 1
@inline xy(i::Int, L::Int) = ((i-1) ÷ L + 1, (i-1) % L + 1)

# symmetric difference of two sets given as integer vectors
function symdiff(a::Vector{Int}, b::Vector{Int})
    seen = Dict{Int,Int}()
    @inbounds for v in a
        seen[v] = get(seen, v, 0) ⊻ 1
    end
    @inbounds for v in b
        seen[v] = get(seen, v, 0) ⊻ 1
    end
    [k for (k,bit) in seen if bit == 1]
end

# build reference "Z" loops as sets of vertex indices
# Cx: winds along x (row at fixed y0); Cy: winds along y (column at fixed x0)
function build_Cx(L::Int; y0::Int=1)
    y = ((y0-1) % L) + 1
    [idx(x, y, L) for x in 1:L]
end
function build_Cy(L::Int; x0::Int=1)
    x = ((x0-1) % L) + 1
    [idx(x, y, L) for y in 1:L]
end

# parity on an indicator set
# @inline function parity_on(set_vec::Vector{Int}, support_vec::Vector{Int})
#     supp = Set(support_vec)
#     c = 0
#     @inbounds for v in set_vec
#         c ⊻= (v in supp) ? 1 : 0
#     end
#     c == 0
# end
@inline function parity_on(S::Vector{Int}, support::Vector{Int})
    supp = BitSet(support)  # fast membership
    p = 0
    @inbounds for v in S
        if v in supp
            p ⊻= 1
        end
    end
    return p == 0
end

# compute plaquette syndrome (Z-checks) for a vertex error set S on an LxL torus
# Each plaquette at lower-left corner (x,y) touches vertices:
# (x,y), (x+1,y), (x,y+1), (x+1,y+1)
function plaquette_syndrome_zero(S::Vector{Int}, L::Int)
    occ = BitSet(S)  # membership set of active vertices
    @inbounds for x in 1:L, y in 1:L
        v1 = idx(x,   y,   L)
        v2 = idx(x+1, y,   L)
        v3 = idx(x,   y+1, L)
        v4 = idx(x+1, y+1, L)
        s = (Int(v1 in occ) + Int(v2 in occ) + Int(v3 in occ) + Int(v4 in occ)) % 2
        if s == 1
            return false
        end
    end
    return true
end

function debug_plaquette_syndrome(S::Vector{Int}, L::Int)
    occ = BitSet(S)
    bad = Tuple{Int,Int}[]
    @inbounds for x in 1:L, y in 1:L
        v1 = idx(x,   y,   L)
        v2 = idx(x+1, y,   L)
        v3 = idx(x,   y+1, L)
        v4 = idx(x+1, y+1, L)
        s = (Int(v1 in occ) + Int(v2 in occ) + Int(v3 in occ) + Int(v4 in occ)) % 2
        if s == 1
            push!(bad, (x,y))
        end
    end
    bad
end


"""
    check_decode(errors::Vector{Int}, errors_true::Vector{Int}, L::Int; x0=1, y0=1)
→ (syndrome_zero::Bool, homology_trivial::Bool)

1) Forms the difference pattern Sdiff = errors ⊕ errors_true (mod-2).
2) Checks that all plaquette Z-checks have even parity on Sdiff.
3) Checks homology triviality by testing parity against two fixed non-contractible Z loops:
   - w_x detected by column Cy(x0)
   - w_y detected by row    Cx(y0)
Returns (syndrome_zero, (w_x==0 && w_y==0)).
"""
# function check_decode(errors::Vector{Int}, errors_true::Vector{Int}, L::Int; x0::Int=1, y0::Int=1)
#     Sdiff = symdiff(errors, errors_true)

#     # (1) plaquette syndrome zero
#     syn0 = plaquette_syndrome_zero(Sdiff, L)

#     # (2) homology triviality via intersection parities
#     Cx = build_Cx(L; y0=y0)  # row loop (winds along x)
#     Cy = build_Cy(L; x0=x0)  # column loop (winds along y)

#     wx_zero = parity_on(Sdiff, Cy)  # true iff no x-winding
#     wy_zero = parity_on(Sdiff, Cx)  # true iff no y-winding

#     return syn0, (wx_zero && wy_zero)
# end
function check_decode(errors::Vector{Int}, errors_true::Vector{Int}, L::Int; x0::Int=1, y0::Int=1)
    Sdiff = symdiff(errors, errors_true)

    # (1) plaquette syndrome zero (already empty-safe)
    syn0 = plaquette_syndrome_zero(Sdiff, L)

    # (2) homology: parities against fixed noncontractible loops
    Cx = build_Cx(L; y0=y0)  # row loop (winds along x)
    Cy = build_Cy(L; x0=x0)  # column loop (winds along y)

    wx_zero = parity_on(Sdiff, Cy)  # x-winding detected by Cy
    wy_zero = parity_on(Sdiff, Cx)  # y-winding detected by Cx

    return syn0, (wx_zero && wy_zero)
end


using Random

# ---------- helpers (same conventions as before) ----------
@inline idx(x::Int, y::Int, L::Int) = (x-1) % L * L + (y-1) % L + 1

function symdiff(a::Vector{Int}, b::Vector{Int})
    seen = Dict{Int,Int}()
    @inbounds for v in a
        seen[v] = get(seen, v, 0) ⊻ 1
    end
    @inbounds for v in b
        seen[v] = get(seen, v, 0) ⊻ 1
    end
    [k for (k,bit) in seen if bit == 1]
end

# four vertices around the plaquette with lower-left corner (x,y)
plaquette_vertices(x::Int, y::Int, L::Int) = [
    idx(x,   y,   L),
    idx(x+1, y,   L),
    idx(x,   y+1, L),
    idx(x+1, y+1, L),
]

# random errors_true: include each vertex with Bernoulli(p)
function random_errors_true(L::Int; p::Float64=0.15, rng=Random.default_rng())
    S = Int[]
    @inbounds for x in 1:L, y in 1:L
        if rand(rng) < p
            push!(S, idx(x,y,L))
        end
    end
    S
end

# pick a random plaquette, XOR it into errors_true to get errors_dec
function add_random_plaquette_mod2(errors_true::Vector{Int}, L::Int; rng=Random.default_rng())
    x = rand(rng, 1:L)
    y = rand(rng, 1:L)
    loop = plaquette_vertices(x,y,L)
    errors_dec = symdiff(errors_true, loop)
    return errors_dec, (x,y), loop
end


# Safer plaquette syndrome: sum Booleans and take mod 2
function plaquette_syndrome_zero(S::Vector{Int}, L::Int)
    occ = BitSet(S)  # membership set of active vertices
    @inbounds for x in 1:L, y in 1:L
        v1 = idx(x,   y,   L)
        v2 = idx(x+1, y,   L)
        v3 = idx(x,   y+1, L)
        v4 = idx(x+1, y+1, L)
        s = (Int(v1 in occ) + Int(v2 in occ) + Int(v3 in occ) + Int(v4 in occ)) % 2
        if s == 1
            return false
        end
    end
    return true
end

# Optional: show which plaquettes are odd
function debug_plaquette_syndrome(S::Vector{Int}, L::Int)
    occ = BitSet(S)
    bad = Tuple{Int,Int}[]
    @inbounds for x in 1:L, y in 1:L
        v1 = idx(x,   y,   L)
        v2 = idx(x+1, y,   L)
        v3 = idx(x,   y+1, L)
        v4 = idx(x+1, y+1, L)
        s = (Int(v1 in occ) + Int(v2 in occ) + Int(v3 in occ) + Int(v4 in occ)) % 2
        if s == 1
            push!(bad, (x,y))
        end
    end
    bad
end

# ------------- tiny sanity check -------------
if abspath(PROGRAM_FILE) == @__FILE__
    L = 5

    # (true, true) — perfect match, no syndrome, trivial homology
    errors_true = Int[]
    errors_dec  = Int[]
    println(check_decode(errors_dec, errors_true, L))  # (true, true)

    # (true, false) — nontrivial loop, no syndrome
    loop_row = build_Cx(L; y0=1)
    errors_bad = vcat(errors_dec, loop_row)
    println(check_decode(errors_bad, errors_true, L))  # (true, false)

    # (false, true) — OPEN pattern with even parity on both reference loops
    # Pick two vertices neither on row y0=1 nor producing odd column parity:
    # e.g., both on column x=1 so column parity is even (2), and neither on row y=1.
    open_trivial = [idx(1,2,L), idx(1,4,L)]
    println(check_decode(open_trivial, errors_true, L))  # (false, true)

    # (false, false) — open string with winding in x-direction
    open_wrap = vcat(open_trivial, loop_row)
    println(check_decode(open_wrap, errors_true, L))  # (false, false)

    rng = MersenneTwister(42)
    L = 7

    # errors_true = random_errors_true(L; p=0.20, rng)
    # errors_dec, (xp,yp), loop = add_random_plaquette_mod2(errors_true, L; rng)

    # println("Random plaquette added at (x,y)=($xp,$yp)")
    # println("sizes: |errors_true|=$(length(errors_true)), |loop|=$(length(loop)), |errors_dec|=$(length(errors_dec))")

    # # Since a plaquette loop is a boundary, the diff is exactly that boundary:
    # # syndrome should stay zero and homology should stay trivial.
    # println("check_decode(errors_dec, errors_true, L) => ", check_decode(errors_dec, errors_true, L))
    # # expected: (true, true)
    errors_true = [rand(1:L^2) for _ in 1:12]              # any 12 vertices
    plaquette = plaquette_vertices(1,5,L)                  # your (x,y) = (1,5)
    errors_dec = symdiff(errors_true, plaquette)

    # sanity: the XOR diff must be exactly that plaquette
    Sdiff = symdiff(errors_dec, errors_true)
    println("Sdiff == plaquette?  ", Set(Sdiff) == Set(plaquette))

    println("check_decode => ", check_decode(errors_dec, errors_true, L))
    println("bad plaquettes: ", debug_plaquette_syndrome(Sdiff, L))
end
