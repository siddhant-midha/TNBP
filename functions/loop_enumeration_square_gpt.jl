# Julia ≥ 1.9

########## Types and basic utilities ##########

const Coord = NTuple{2,Int}                       # (x,y)
const Edge  = NTuple{2,Coord}                     # ((x1,y1),(x2,y2))

# Lexicographic order for coordinates
@inline leq_coord(a::Coord, b::Coord) = (a[1] < b[1]) || (a[1] == b[1] && a[2] <= b[2])

# Normalized, ordered edge
@inline norm_edge(u::Coord, v::Coord)::Edge = leq_coord(u,v) ? (u,v) : (v,u)

# Four-neighborhood in the (unwrapped) square lattice
@inline neighbors(v::Coord) = ((v[1]+1, v[2]), (v[1]-1, v[2]), (v[1], v[2]+1), (v[1], v[2]-1))

# Turn an edge-set into a sortable string key
function edge_list_key(E::Set{Edge})
    es = collect(E)
    sort!(es, by = e -> (e[1][1], e[1][2], e[2][1], e[2][2]))
    io = IOBuffer()
    for e in es
        print(io, e[1][1], ',', e[1][2], ';', e[2][1], ',', e[2][2], '|')
    end
    return String(take!(io))
end

########## D4 symmetries about the root (0,0) ##########

# 8 lattice symmetries that fix the origin
@inline T_id(p::Coord)      = ( p[1],  p[2])
@inline T_r90(p::Coord)     = (-p[2],  p[1])
@inline T_r180(p::Coord)    = (-p[1], -p[2])
@inline T_r270(p::Coord)    = ( p[2], -p[1])
@inline T_refx(p::Coord)    = (-p[1],  p[2])    # reflect across y-axis
@inline T_refy(p::Coord)    = ( p[1], -p[2])    # reflect across x-axis
@inline T_refdiag(p::Coord) = ( p[2],  p[1])    # reflect across y=x
@inline T_refanti(p::Coord) = (-p[2], -p[1])    # reflect across y=-x

const D4 = (T_id, T_r90, T_r180, T_r270, T_refx, T_refy, T_refdiag, T_refanti)

# Canonical key of an edge-set under D4 (root fixed at origin)
function canonical_key(E::Set{Edge})
    best::Union{Nothing,String} = nothing
    for T in D4
        # transform and renormalize each edge
        E2 = Set{Edge}()
        for (u,v) in E
            push!(E2, norm_edge(T(u), T(v)))
        end
        k = edge_list_key(E2)
        if best === nothing || k < best
            best = k
        end
    end
    return best::String
end

# Canonicalized edge-set itself (useful for output)
function canonical_edges(E::Set{Edge})
    best_key::Union{Nothing,String} = nothing
    best_E::Union{Nothing,Set{Edge}} = nothing
    for T in D4
        E2 = Set{Edge}()
        for (u,v) in E
            push!(E2, norm_edge(T(u), T(v)))
        end
        k = edge_list_key(E2)
        if best_key === nothing || k < best_key
            best_key, best_E = k, E2
        end
    end
    return best_E::Set{Edge}
end

########## State + helpers ##########

mutable struct State
    V::Set{Coord}               # vertices present
    E::Set{Edge}                # edges present
    deg::Dict{Coord,Int}        # current degrees within the subgraph
end

# Compute the frontier: all lattice edges with ≥1 endpoint in V that are not yet in E
function frontier_edges(V::Set{Coord}, E::Set{Edge})
    F = Set{Edge}()
    for v in V
        for u in neighbors(v)
            e = norm_edge(v,u)
            if !(e in E)
                push!(F, e)
            end
        end
    end
    return F
end

# Add an edge and update degrees/vertex set
function add_edge(s::State, e::Edge)
    (u,v) = e
    V2 = copy(s.V)
    E2 = copy(s.E)
    deg2 = copy(s.deg)
    push!(E2, e)
    if !(u in V2); push!(V2,u); deg2[u] = 0; end
    if !(v in V2); push!(V2,v); deg2[v] = 0; end
    deg2[u] = get(deg2,u,0) + 1
    deg2[v] = get(deg2,v,0) + 1
    return State(V2, E2, deg2)
end

# Sum of degree deficits (how many half-edges still needed to reach degree ≥ 2 everywhere)
@inline total_deficit(deg::Dict{Coord,Int}) = sum(d < 2 ? 2 - d : 0 for d in values(deg))

########## BFS enumeration with early pruning and canonicalization ##########

"""
    enumerate_loops(max_edges; root=(0,0), min_edges=2)

Enumerate connected “loops” (each vertex has degree ≥ 2) on the 2D square lattice,
supported on `root`, by growing edges with BFS. Uses:

- Early pruning: if `remaining_edges < ceil(deficit/2)`, prune.
- Canonicalization under the D4 symmetries about the root.

Returns a `Dict{Int, Vector{Vector{Edge}}}` mapping edge-count → list of canonicalized loops.
Each loop is a sorted `Vector{Edge}` with edges as `((x1,y1),(x2,y2))`.
"""
function enumerate_loops(max_edges::Int; root::Coord=(0,0), min_edges::Int=2)
    # Initial state: only the root present, no edges, degree 0
    init = State(Set([root]), Set{Edge}(), Dict(root => 0))

    # BFS queue and head pointer
    queue = Vector{State}([init])
    head = 1

    # Seen canonical keys per edge-count (prevents re-enqueue duplicates)
    seen = Dict{Int, Set{String}}()
    seen[0] = Set([canonical_key(init.E)])

    # Collected results by size
    results = Dict{Int, Vector{State}}()

    while head <= length(queue)
        s = queue[head]; head += 1
        m = length(s.E)

        # Canonical check (expand only the canonical representatives of partial states)
        if edge_list_key(s.E) != canonical_key(s.E)
            continue
        end

        # Early pruning: each new edge can reduce deficit by at most 2
        deficit = total_deficit(s.deg)
        remaining = max_edges - m
        if remaining < cld(deficit, 2)   # ceil(deficit/2)
            continue
        end

        # If this state is already a valid loop, record it (we may still expand further)
        if deficit == 0 && m >= min_edges
            push!(get!(results, m, Vector{State}()), s)
        end

        # Do not expand beyond max_edges
        if m == max_edges
            continue
        end

        # Expand by adding each frontier edge once
        for e in frontier_edges(s.V, s.E)
            s2 = add_edge(s, e)
            m2 = m + 1
            cank = canonical_key(s2.E)
            bucket = get!(seen, m2, Set{String}())
            if !(cank in bucket)
                push!(bucket, cank)
                push!(queue, s2)
            end
        end
    end

    # Canonicalize outputs (edge-sets) and deduplicate within each size
    out = Dict{Int, Vector{Vector{Edge}}}()
    for (k, vec) in results
        uniq = Set{String}()
        acc  = Vector{Vector{Edge}}()
        for s in vec
            Ecan = canonical_edges(s.E)
            key  = edge_list_key(Ecan)
            if !(key in uniq)
                push!(uniq, key)
                es = collect(Ecan)
                sort!(es, by = e -> (e[1][1], e[1][2], e[2][1], e[2][2]))
                push!(acc, es)
            end
        end
        out[k] = acc
    end
    return out
end

########## Example usage ##########
# Enumerate all loops with up to 8 edges supported at (0,0):
# loops = enumerate_loops(8; root=(0,0), min_edges=4)
# println("sizes: ", sort!(collect(keys(loops))))
# println("count at 6 edges: ", length(get(loops, 6, [])))
