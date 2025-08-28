include("tc_mld.jl")


function LER(L, p; nsamples = 100)
    n = 2 * L^2 
    pcmat = toric_code_X_parity_matrix(L)
    logicals = toric_code_logical_operators(L)
    logical_error_list = []
    for nsamp = 1:nsamples
        errors_true = [rand() < p ? 1 : 0 for _ in 1:n]
        # Compute syndrome
        syndrome = pcmat * errors_true .% 2
        # Get a consistent error configuration
        matched_errors = simple_match_syndrome(syndrome, L)

        # Get logical data (which logical operators are applied)
        logical_data = (logicals * matched_errors) .% 2

        # Create trivial error chain (remove logical component from matched_errors)
        trivial_error_chain = (matched_errors .+ (logical_data[1] .* logicals[1,:]) .+ (logical_data[2] .* logicals[2,:])) .% 2

        # Generate all four logical equivalence classes
        error_chains = [
            trivial_error_chain,                                                  # 00
            (trivial_error_chain .+ logicals[1,:]) .% 2,                          # 10  
            (trivial_error_chain .+ logicals[2,:]) .% 2,                          # 01
            (trivial_error_chain .+ logicals[1,:] .+ logicals[2,:]) .% 2          # 11
        ]
        β = β_(p)
        free_energies = [BP_free_energy(error_chain, L, β) for error_chain in error_chains]
        logical_error = Int(sum((logicals * error_chains[argmax(free_energies)]) .% 2) > 0)
        push!(logical_error_list, logical_error)
    end 
    return mean(logical_error_list)
end 

L = 5
p_ar = 0:0.01:0.1  
nsamples = 20
ler_ar = [LER(L, p; nsamples = nsamples) for p in p_ar]

pl = plot(p_ar, ler_ar, marker=:square)

display(pl)

readline()