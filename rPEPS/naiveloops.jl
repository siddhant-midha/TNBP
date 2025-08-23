include("randompeps.jl")


ti = true 
orthog = false 
annealing = 0.75 
maxiter = 1000 


## i will do it for a fixed (N,w,\eta) and one sample 

N = 5
T = N
w = 4 
η = 0.1
nsamples = 1

cluster_data = load_latest_cluster_file(N,w)
loop_objects = cluster_data.all_loops  
all_loops = [loop_object.edges for loop_object in loop_objects]


tensors, peps = peps_controllable(N, T; η=η, ti=ti, orthog=orthog)
exact_FE = log(real(contract_peps_no_phys(peps; cutoff=1E-12, maxdim=2^N)))

adj_mat, edges, links = BP.get_adj_mat(tensors)
messages = BP.get_messages(tensors, edges, links) 
messages = BP.message_passing(tensors, messages, edges, adj_mat; α=annealing, max_iters=maxiter, diagnose=false, normalise=true)
Z_list = BP.get_fixed_point_list(tensors, messages, adj_mat)
bp_FE = sum(log.(real.(Z_list)))
tensors = BP.normalize_tensors(tensors, Z_list)

global loopcorrx = 1.

for loop in all_loops
    if all([e in edges for e in loop])
        loop_contri = scalar(BP.loop_contribution(loop, messages, tensors, edges, links, adj_mat))
        global loopcorrx += loop_contri
    end
end


println("BP error: ", abs(bp_FE - exact_FE))
println("BP + loop correction error: ", abs(bp_FE + log(loopcorrx) - exact_FE))