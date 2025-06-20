include("dependencies.jl")
include("brute_force.jl")
include("enumerate_loops.jl")
include("boundary_evolution.jl")

# test loop
let
    run(`clear`)
    loop_l = loop10th(5,5,6,7)
    # show number of loops
    @show length(loop_l)
    for loop in loop_l
        # test no open edges and at most degree four
        all_vertices_have_multiple_edges = all(v -> (degree(loop, v) != 1 && degree(loop, v) <= 4), vertices(loop))
        if !all_vertices_have_multiple_edges
            @show loop
        end
        # test number of edges
        if ne(loop) != 10
            @show loop
        end
    end
end

# test loop enumeration
let 
    run(`clear`)
    g = square_lattice_periodic(7,7)
    m = 4
    loops = loops_on_vertex_edge_growth(g, 1, m)
    
    println("Found $(length(loops)) unique loops (≤$m edges) supported on vertex 1")
    for (k,ℓ) in enumerate(loops)
        # println("  loop $k has $(ne(ℓ)) edges")

        if ne(ℓ) > m
            @show ℓ
        end
    end

    loop = loops[1]
    @show loop
    # for edge in edges(loop)
    #     println("Edge: ", edge)
    # end
    # Show the list of edges and the number of vertices
    # println("Edges: ", edges(loop))
    # println("Number of vertices: ", length(vertices(loop)))
end

let 
    run(`clear`)
    g = square_lattice_periodic(6, 6)
    m = 6
    loops = all_loops(g, m)
    println("unique connected loops (≤10 edges) in 6×6 torus: $(length(loops))")
    
end

# test boundary evolution
let
    run(`clear`)
    
    # Define the Ising model parameters
    Ny = 10 # Height
    Nx = 10 # Width
    beta = 0.4 # Inverse temperature (critical beta is ~0.4406)
    hz = 0.0   # External magnetic field

    # Generate the PEPS for the Ising model partition function
    peps = generate_ising_peps(Ny, Nx, beta, hz)

    # Perform the approximate contraction
    partition_function = contract_peps_no_phys(peps; cutoff=1E-8, maxdim=40)

    # println("\nApproximate Partition Function Z: ", partition_function)
    println("\nFree Energy Density: ", log(partition_function) / (Ny * Nx))
    
    # For a 4x4 lattice at beta=0.4, hz=0, the exact result is ~2.367e5
    # The approximate result should be close depending on `maxdim`.
end