include("dependencies.jl")
include("brute_force.jl")

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