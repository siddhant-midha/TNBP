const Edge = Tuple{Int,Int}
edge(u, v) = u < v ? (u, v) : (v, u)
make_graph(nv, E::Set{Edge}) = (g = SimpleGraph(nv); foreach(e->add_edge!(g,e...), E); g)
verts_of_edges(E) = Set(v for (u,v) in E for v in (u,v))


# --- NEW: keep vertex identities; two loops are the same
#          iff they have exactly the same edge set ------------------
function edgekey(E::Set{Edge})
    join(sort(collect(E)), ";")      # canonical string for THIS labelled graph
end

# ------------------------------------------------------------------
#  Weisfeiler-Lehman vertex-colour refinement → canonical hash
# ------------------------------------------------------------------
"""
    wl_hash(edge_set::Set{Tuple{Int,Int}}, n::Int; rounds::Int = 3)

Return a WL hash string of an *edge-induced* simple graph on `n` vertices
(labelled `1:n`) whose edges are `edge_set`.  
Isomorphic graphs produce the same hash with very high probability.
"""
function wl_hash(edge_set::Set{Edge}, n::Int; rounds::Int=3)
    # adjacency lists
    adj = [Int[] for _ in 1:n]
    for (u, v) in edge_set
        push!(adj[u], v)
        push!(adj[v], u)
    end

    # initial colours = degrees
    colour = [length(adj[v]) for v in 1:n]

    for _ in 1:rounds
        # new colour = hash(old_colour, multiset(neighbour colours))
        new_raw = Vector{UInt64}(undef, n)
        for v in 1:n
            neigh = sort(colour[adj[v]])
            new_raw[v] = hash((colour[v], neigh))
        end
        # re-index to dense integers for stability
        remap = Dict{UInt64,Int}()
        next_id = 0
        for h in new_raw
            haskey(remap, h) || (next_id += 1; remap[h] = next_id)
        end
        colour = [remap[h] for h in new_raw]
    end
    return join(sort(colour), ',')       # order-independent signature
end

# ------------------------------------------------------------------
#  Edge-growth search for loops
# ------------------------------------------------------------------
all_deg_gt1(deg) = all(>(1), values(deg))

function subgraph_from_edges(n::Int, edge_set::Set{Edge})
    sg = SimpleGraph(n)
    for (u, v) in edge_set
        add_edge!(sg, u, v)
    end
    return sg
end

"""
    loops_on_vertex_edge_growth(g::SimpleGraph, root::Int, m::Int)

Return all *unique* loops (edge-induced, connected, every vertex deg>1)
containing `root` and using ≤ `m` edges.  
Each result is a `SimpleGraph` with **the same vertex count as `g`**; vertices
not in the loop remain isolated.
"""
function loops_on_vertex_edge_growth(g::SimpleGraph, root::Int, m::Int)
    n = nv(g)
    loops = SimpleGraph[]
    seen = Set{String}()

    function dfs(edges::Set{Edge},
        deg::Dict{Int,Int},
        verts::Set{Int},
        remaining::Int)

        if length(edges) ≥ 3 && all_deg_gt1(deg)
            # h = wl_hash(edges, n)
            h = edgekey(edges)
            if h ∉ seen
                push!(seen, h)
                push!(loops, subgraph_from_edges(n, edges))
            end
            # push!(seen, h)
            # push!(loops, subgraph_from_edges(n, edges))
        end
        remaining == 0 && return

        # candidate frontier edges
        front = Set{Edge}()
        for u in verts, v in neighbors(g, u)
            e = edge(u, v)
            if e ∉ edges && length(edges) + 1 ≤ m
                push!(front, e)
            end
        end

        for (u, v) in front
            new_edges = copy(edges)
            push!(new_edges, (u, v))
            new_verts = verts ∪ Set([u, v])

            new_deg = copy(deg)
            new_deg[u] = get(new_deg, u, 0) + 1
            new_deg[v] = get(new_deg, v, 0) + 1

            # pruning: each degree-1 vertex needs at least one more edge
            need = count(==(1), values(new_deg))
            remaining - 1 < need && continue

            dfs(new_edges, new_deg, new_verts, remaining - 1)
        end
    end

    dfs(Set{Edge}(), Dict(root => 0), Set([root]), m + 1)
    return loops
end

# ------------------------------------------------------------------
# grow edge-sets from one root
# ------------------------------------------------------------------
function grow_from_root!(g::SimpleGraph, root::Int, m::Int,
    seen::Set{String}, loops::Vector{SimpleGraph})

    n = nv(g)

    function dfs(edges::Set{Edge}, deg::Dict{Int,Int},
        verts::Set{Int}, remaining::Int)

        if length(edges) ≥ 3 && all_deg_gt1(deg)
            sig = edgekey(edges)
            if sig ∉ seen                     # <- global deduplication
                push!(seen, sig)
                push!(loops, make_graph(n, edges))
            end
        end
        remaining == 0 && return

        # frontier = unused edge touching current vertices
        front = Set{Edge}()
        for u in verts, v in neighbors(g, u)
            e = edge(u, v)
            if e ∉ edges && length(edges) + 1 ≤ m
                push!(front, e)
            end
        end

        for (u, v) in front
            new_edges = copy(edges)
            push!(new_edges, (u, v))
            new_verts = verts ∪ (u == v ? Set([u]) : Set([u, v]))
            new_deg = copy(deg)
            new_deg[u] = get(new_deg, u, 0) + 1
            new_deg[v] = get(new_deg, v, 0) + 1

            # prune if remaining edges can’t lift every deg-1 vertex
            need = count(==(1), values(new_deg))
            remaining - 1 < need && continue

            dfs(new_edges, new_deg, new_verts, remaining - 1)
        end
    end

    dfs(Set{Edge}(), Dict(root => 0), Set([root]), m+1)
end

# ------------------------------------------------------------------
# enumerate _all_ connected loops of weight ≤ m, removing duplicates
# ------------------------------------------------------------------
function all_loops(g::SimpleGraph, m::Int)
    n         = nv(g)
    seen      = Set{String}()          # ← GLOBAL for the whole sweep
    unique_loops = SimpleVector()      # will hold each loop once

    function dfs(E::Set{Tuple{Int,Int}},
                 deg::Dict{Int,Int},
                 V::Set{Int},
                 remaining::Int)

        # record a loop if new
        if length(E) ≥ 3 && all(>(1), values(deg))
            sig = edgekey(E)
            if sig ∉ seen                      # duplicates filtered here
                push!(seen, sig)
                push!(unique_loops,
                      ( graph  = (sg = SimpleGraph(n);
                                  foreach(e -> add_edge!(sg, e...), E); sg),
                        edges  = copy(E),
                        verts  = copy(V),
                        weight = length(E) ))
            end
        end
        remaining == 0 && return

        # frontier = unused edge touching the current vertex set
        front = Set{Tuple{Int,Int}}()
        for u in V, v in neighbors(g, u)
            e = u < v ? (u,v) : (v,u)
            length(E) + 1 ≤ m && e ∉ E && push!(front, e)
        end

        for (u,v) in front
            newE   = copy(E);  push!(newE, (u,v))
            newV   = V ∪ Set([u,v])
            newdeg = copy(deg);
            newdeg[u] = get(newdeg, u, 0) + 1
            newdeg[v] = get(newdeg, v, 0) + 1
            need = count(==(1), values(newdeg))
            remaining - 1 < need && continue
            dfs(newE, newdeg, newV, remaining - 1)
        end
    end

    # -------- sweep every vertex as a seed, reuse ONE `seen` set ----
    for v in vertices(g)
        dfs(Set{Tuple{Int,Int}}(), Dict(v => 0), Set([v]), m)
    end
    return collect(unique_loops)   # each loop appears exactly once
end

# ------------------------------------------------------------------
# 2.  interaction graph of loops
# ------------------------------------------------------------------
function interaction_graph(loop_verts::Vector{Set{Int}})
    L = length(loop_verts)
    G = SimpleGraph(L)
    for i in 1:L-1, j in i+1:L
        !isempty(loop_verts[i] ∩ loop_verts[j]) && add_edge!(G,i,j)
    end
    G
end

# ------------------------------------------------------------------
#  connectivity test for a multi-loop
# ------------------------------------------------------------------
function connected_union(nv, edge_sets::Vector{Set{Edge}})
    U = SimpleGraph(nv)
    for E in edge_sets
        foreach(e->add_edge!(U,e...), E)
    end
    Vsub = [v for v in 1:nv if degree(U,v) > 0]
    length(Vsub) ≤ 1 && return true
    return is_connected(induced_subgraph(U, Vsub))
end

# ------------------------------------------------------------------
# 3.  enumerate connected multi-loops by growing a connected node set
# ------------------------------------------------------------------
"""
    multiloops_supported_on_site(g, root, m)

Return every **connected multi-loop** P with total weight ≤ m containing
≥ 1 loop that passes through `root`.  Each P is
`Vector{Tuple{SimpleGraph,Int}}`.
"""
function multiloops_supported_on_site(g::SimpleGraph, root::Int, m::Int)

    # ----- single loops & metadata
    loops = all_loops(g, m)
    L     = length(loops)
    weights  = [ℓ.weight for ℓ in loops]
    verts    = [ℓ.verts  for ℓ in loops]
    root_idx = [i for i in 1:L if root in verts[i]]
    isempty(root_idx) && error("no loop contains root")

    # ----- interaction graph
    IG = interaction_graph(verts)
    neigh = [neighbors(IG,i) for i in 1:L]

    # ----- DFS over connected node sets + multiplicities
    results, seenP = Vector{Vector{Tuple{SimpleGraph,Int}}}(), Set{String}()

    function sig_multiset(idx::Vector{Int}, mult::Vector{Int})
        join(sort(["$(idx[i])^$(mult[i])" for i in eachindex(idx)]), "|")
    end

    function dfs(nodes::Vector{Int}, mult::Vector{Int},
                 weight::Int, frontier::Set{Int})

        # record if root included and under budget
        if weight ≤ m && !isempty(nodes) && any(i->i in root_idx, nodes)
            key = sig_multiset(nodes, mult)
            if key ∉ seenP
                push!(seenP, key)
                push!(results, [ (loops[nodes[i]].graph, mult[i])
                                  for i in eachindex(nodes) ])
            end
        end
        weight == m && return   # cannot add more

        # ---- (A) increase multiplicity of an existing node
        for i in eachindex(nodes)
            w = weights[nodes[i]]
            weight + w > m && continue
            mult[i] += 1
            dfs(nodes, mult, weight+w, frontier, )
            mult[i] -= 1
        end

        # ---- (B) add a new frontier node to keep connectivity
        for v in copy(frontier)
            w = weights[v]
            weight + w > m && continue
            push!(nodes, v);  push!(mult, 1)
            new_front = (frontier ∪ neigh[v]) ∖ Set(nodes)
            dfs(nodes, mult, weight+w, new_front)
            pop!(nodes); pop!(mult)
        end
    end

    # start: pick each root-containing loop as seed with mult = 1..max
    for r in root_idx
        maxk = m ÷ weights[r]
        for k in 1:maxk
            dfs([r], [k], k*weights[r], Set(neigh[r]))
        end
    end
    results
end


# ------------------------------------------------------------------
#  Tiny demo on a 4×4 square lattice
# ------------------------------------------------------------------
function square_lattice(n, m)
    g = SimpleGraph(n * m)
    for i in 1:n, j in 1:m
        v = (i - 1) * m + j
        j < m && add_edge!(g, v, v + 1)
        i < n && add_edge!(g, v, v + m)
    end
    g
end


"""
    square_lattice_periodic(nx, ny) -> SimpleGraph

`nx × ny` grid with **periodic boundaries** (torus).  
Vertices are labelled 1 … nx·ny in row-major order.
"""
function square_lattice_periodic(nx::Int, ny::Int)
    g = SimpleGraph(nx * ny)

    # helper: convert 2-D coords to 1-D label (1-based)
    index(i, j) = (i % nx) * ny + (j % ny) + 1   # (% is mod, 0-based then +1)

    for i in 0:nx-1, j in 0:ny-1
        v = index(i, j)
        vr = index(i, j + 1)   # right neighbour  (wraps on j)
        vd = index(i + 1, j)   # down  neighbour  (wraps on i)
        add_edge!(g, v, vr)
        add_edge!(g, v, vd)
    end
    return g
end