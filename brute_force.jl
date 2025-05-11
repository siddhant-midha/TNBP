function coordToIdx(i,j,dx,dy)
    (i - 1) * dy + j
end

# return a list of connected loops with origin (origin_x, origin_y) and has weight 4 
# (origin_x, origin_y): origin of the loop, defined as the lower-left corner
# (dx, dy): size of the PEPS
# return: an array of Graph objects
function loop4th(origin_x, origin_y, dx,dy)
    loop_l = []
    
    loop = Graph(dx*dy)

    start_x = origin_x
    start_y = origin_y
    end_x = mod1(origin_x+1,dx)
    end_y = origin_y
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = origin_y
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = origin_x
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = origin_x
    start_y = mod1(origin_y+1,dy)
    end_x = origin_x
    end_y = origin_y
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)
    return loop_l
end

# return a list of connected loops with origin (origin_x, origin_y) and has weight 6
# (origin_x, origin_y): origin of the loop, defined as the lower-left corner
# (dx, dy): size of the PEPS
# return: an array of Graph objects
function loop6th(origin_x, origin_y, dx,dy)
    loop_l = []
    
    # horizontal
    loop = Graph(dx*dy)

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # vertical
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    return loop_l
end

# return a list of connected loops with origin (origin_x, origin_y) and has weight 7
# (origin_x, origin_y): origin of the loop, defined as the lower-left corner
# (dx, dy): size of the PEPS
# return: an array of Graph objects
function loop7th(origin_x, origin_y, dx,dy)
    loop_l = []
    
    # horizontal
    loop = Graph(dx*dy)

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # vertical
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    return loop_l
end

# return a list of connected loops with origin (origin_x, origin_y) and has weight 8
# (origin_x, origin_y): origin of the loop, defined as the lower-left corner
# (dx, dy): size of the PEPS
# return: an array of Graph objects
function loop8th(origin_x, origin_y, dx,dy)
    loop_l = []
    
    # strip horizontal
    loop = Graph(dx*dy)

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # strip vertical
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair SW
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair SE
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair NE
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair NW
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    return loop_l
end

# return a list of connected loops with origin (origin_x, origin_y) and has weight 9
# (origin_x, origin_y): origin of the loop, defined as the lower-left corner
# (dx, dy): size of the PEPS
# return: an array of Graph objects
function loop9th(origin_x, origin_y, dx,dy)
    loop_l = []

    # strip W
    loop = Graph(dx*dy)

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # strip E
    loop = Graph(dx*dy)

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # strip N
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # strip S
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair SW clockwise
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair SE clockwise
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair NE clockwise
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair NW clockwise
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair SW counter-clockwise
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair SE counter-clockwise
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair NE counter-clockwise
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair NW counter-clockwise
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    return loop_l
end

# return a list of connected loops with origin (origin_x, origin_y) and has weight 10
# (origin_x, origin_y): origin of the loop, defined as the lower-left corner
# (dx, dy): size of the PEPS
# return: an array of Graph objects
function loop10th(origin_x, origin_y, dx,dy)
    loop_l = []

    # strip horizontal
    loop = Graph(dx*dy)

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # strip vertical
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair SW
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair SE
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair NE
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # chair NW
    loop = Graph(dx*dy)
    
    start_x = mod1(origin_x,dx)
    start_y = mod1(origin_y,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # cube SW
    loop = Graph(dx*dy)

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # cube SE
    loop = Graph(dx*dy)

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # cube NE
    loop = Graph(dx*dy)

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # cube NW
    loop = Graph(dx*dy)

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # long strip horizontal
    loop = Graph(dx*dy)

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+4,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+4,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+4,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+4,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+3,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+3,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+2,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+2,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    # long strip vertical
    loop = Graph(dx*dy)

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+0,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x+1,dx)
    end_y = mod1(origin_y+4,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+1,dx)
    start_y = mod1(origin_y+4,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+4,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+4,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+3,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+3,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+2,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+2,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+1,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    start_x = mod1(origin_x+0,dx)
    start_y = mod1(origin_y+1,dy)
    end_x = mod1(origin_x+0,dx)
    end_y = mod1(origin_y+0,dy)
    add_edge!(loop,coordToIdx(start_x,start_y,dx,dy),coordToIdx(end_x,end_y,dx,dy))

    push!(loop_l,loop)

    return loop_l
end

# return the 4th order correction of Z
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function ZCorrection4th(loopValue,dx,dy)
    dZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop4th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dZ += loopValue(loop)
        end
    end

    return dZ
end

# return the 6th order correction of Z
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function ZCorrection6th(dx,dy)
    dZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop6th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dZ += loopValue(loop)
        end
    end

    return dZ
end

# return the 7th order correction of Z
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function ZCorrection7th(dx,dy)
    dZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop7th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dZ += loopValue(loop)
        end
    end

    return dZ
end

# return the 8th order correction of Z
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function ZCorrection8th(dx,dy)
    dZ = 0

    # one loop
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop8th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dZ += loopValue(loop)
        end
    end

    # two loops 4+4
    for origin_y1 = 1:dy, origin_x1 = 1:dx
        for origin_y2 = origin_y1:dy, origin_x1 = origin_x1:dx
            loop1_l = loop4th(origin_x1,origin_y1,dx,dy)
            loop2_l = loop4th(origin_x2,origin_y2,dx,dy)
            for loop1 in loop1_l, loop2 in loop2_l
                # check disconnected
                deg2_loop1 = Set([v for v in vertices(loop1) if degree(loop1, v) ≥ 2])
                deg2_loop2 = Set([v for v in vertices(loop2) if degree(loop2, v) ≥ 2])
                has_common_vertex = !isempty(intersect(deg2_loop1, deg2_loop2))
                if has_common_vertex
                    continue
                end

                dZ += loopValue(loop1)*loopValue(loop2)
            end
        end
    end
    return dZ
end

# return the 9th order correction of Z
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function ZCorrection9th(dx,dy)
    dZ = 0

    # one loop
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop9th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dZ += loopValue(loop)
        end
    end
    return dZ
end

# return the 10th order correction of Z
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function ZCorrection10th(dx,dy)
    dZ = 0

    # one loop
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop10th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dZ += loopValue(loop)
        end
    end

    # two loops 4+6
    for origin_y1 = 1:dy, origin_x1 = 1:dx
        for origin_y2 = origin_y1:dy, origin_x1 = origin_x1:dx
            loop1_l = loop4th(origin_x1,origin_y1,dx,dy)
            loop2_l = loop6th(origin_x2,origin_y2,dx,dy)
            for loop1 in loop1_l, loop2 in loop2_l
                # check disconnected
                deg2_loop1 = Set([v for v in vertices(loop1) if degree(loop1, v) ≥ 2])
                deg2_loop2 = Set([v for v in vertices(loop2) if degree(loop2, v) ≥ 2])
                has_common_vertex = !isempty(intersect(deg2_loop1, deg2_loop2))
                if has_common_vertex
                    continue
                end

                dZ += loopValue(loop1)*loopValue(loop2)
            end
        end
    end
    return dZ
end

# return the 4th order correction of logZ
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function logZCorrection4th(loopValue,dx,dy)
    dlogZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop4th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dlogZ += loopValue(loop)
        end
    end

    return dlogZ
end

# return the 6th order correction of logZ
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function logZCorrection6th(loopValue,dx,dy)
    dlogZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop6th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dlogZ += loopValue(loop)
        end
    end

    return dlogZ
end

# return the 7th order correction of logZ
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function logZCorrection7th(loopValue,dx,dy)
    dlogZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop7th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dlogZ += loopValue(loop)
        end
    end

    return dlogZ
end

# return the 8th order correction of logZ
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function logZCorrection8th(loopValue,dx,dy)
    dlogZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop8th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dlogZ += loopValue(loop)
        end
    end

    # two loops 4+4
    for origin_y1 = 1:dy, origin_x1 = 1:dx
        for origin_y2 = origin_y1:dy, origin_x1 = origin_x1:dx
            loop1_l = loop4th(origin_x1,origin_y1,dx,dy)
            loop2_l = loop4th(origin_x2,origin_y2,dx,dy)
            for loop1 in loop1_l, loop2 in loop2_l
                # check connected
                deg2_loop1 = Set([v for v in vertices(loop1) if degree(loop1, v) ≥ 2])
                deg2_loop2 = Set([v for v in vertices(loop2) if degree(loop2, v) ≥ 2])
                has_common_vertex = !isempty(intersect(deg2_loop1, deg2_loop2))
                if !has_common_vertex
                    continue
                end

                dlogZ -= loopValue(loop1)*loopValue(loop2)/2
            end
        end
    end

    return dlogZ
end

# return the 9th order correction of logZ
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function logZCorrection9th(loopValue,dx,dy)
    dlogZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop9th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dlogZ += loopValue(loop)
        end
    end

    return dlogZ
end

# return the 10th order correction of logZ
# this assumes BP fixed point is normalized to 1
# loopValue(loop): a function that inputs a loop and outputs its loop value
function logZCorrection10th(loopValue,dx,dy)
    dlogZ = 0
    for origin_y = 1:dy, origin_x = 1:dx
        loop_l = loop10th(origin_x,origin_y,dx,dy)
        for loop in loop_l
            dlogZ += loopValue(loop)
        end
    end

    # two loops 4+6
    for origin_y1 = 1:dy, origin_x1 = 1:dx
        for origin_y2 = origin_y1:dy, origin_x1 = origin_x1:dx
            loop1_l = loop4th(origin_x1,origin_y1,dx,dy)
            loop2_l = loop6th(origin_x2,origin_y2,dx,dy)
            for loop1 in loop1_l, loop2 in loop2_l
                # check disconnected
                deg2_loop1 = Set([v for v in vertices(loop1) if degree(loop1, v) ≥ 2])
                deg2_loop2 = Set([v for v in vertices(loop2) if degree(loop2, v) ≥ 2])
                has_common_vertex = !isempty(intersect(deg2_loop1, deg2_loop2))
                if !has_common_vertex
                    continue
                end

                dZ -= loopValue(loop1)*loopValue(loop2)/2
            end
        end
    end
    return dlogZ
end
