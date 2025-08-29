include("randompeps.jl")
using Plots 

function test_2by2_rPEPS(η, type, annealing, maxiter)
    N = 2 
    T = N
    tensors, peps = peps_controllable(N, T; η=η, ti=true, type=type)
    exact_PF = contract_peps_no_phys(peps; cutoff=1E-14, maxdim=2^N)

    # BP calculation - ensure complex arithmetic throughout
    adj_mat, edges, links = BP.get_adj_mat(tensors)
    messages = BP.get_messages(tensors, edges, links; random_part=0.01)
    messages,arr = BP.message_passing(tensors, messages, edges, adj_mat; α=annealing, max_iters=maxiter, diagnose=true)
    p = plot(arr)
    Z_l = BP.get_fixed_point_list(tensors, messages, adj_mat)
    # Ensure Z_l entries are complex before taking log
    Z_l_complex = Complex.(Z_l)
    bp_PF = prod((Z_l_complex)) 

    T_normalized = BP.normalize_tensors(tensors, Z_l)

    corrected_PF = bp_PF * (1 + scalar(BP.loop_contribution(edges, messages, T_normalized, edges, links, adj_mat)))

    println("Type of rPEPS...", type)
    println("Checking self consistency...", BP.check_self_consistency(tensors, messages, adj_mat))
    println()
    println("Exact PF...", exact_PF)
    println()
    println("BP corrected PF...", corrected_PF)
    println()
    println(messages[1,2])
    println()

    display(p)
    readline()
end 

## 2 by 2, BP + loop should give exact answer if converged. Does not converge for complex type even for \eta = 0.


η = 0.0
annealing = 0.5
maxiter = 5000
type = "complex"

test_2by2_rPEPS(η, type, annealing, maxiter)