using ITensors, ITensorMPS, Graphs, LinearAlgebra, Statistics
include("../functions/ClusterEnumeration.jl")
include("../functions/ctmrg.jl")

function load_latest_cluster_file(N, w; bc = "periodic")
    save_dir = "saved_clusters"
    if !isdir(save_dir)
        error("No saved_clusters directory found!")
    end
    files = readdir(save_dir)
    matching_files = filter(f -> endswith(f, ".jld2") && contains(f, "L$N") && contains(f, "w$w") && contains(f, "$bc"), files)
    if isempty(matching_files)
        error("No matching cluster file found for N=$N, w=$w in $save_dir!")
    end
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading cluster data from: $(latest_file)")
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end


function make_pair_combiner(ket::Index, bra::Index, combnd::Index)
    dket, dbra = dim(ket), dim(bra)
    @assert dim(combnd) == dket * dbra "Dimension mismatch: combnd must be dim(ket)*dim(bra)"
    C = ITensor(ket, bra, combnd)
    for n in 1:dket, m in 1:dbra
        C[ket=>n, bra=>m, combnd=>(m + (n-1)*dbra)] = 1.0
    end
    return C
end

function iσy(idx1, idx2)
    ## makes singlets
    tens = ITensor(idx1,idx2)
    tens[1,2] = -1 
    tens[2,1] = 1 
    return tens 
end 

function aklt_square_tensor(;tag_prefix="")
    ## bond index mapping (julia -> actual), 1 -> 0, 2 -> 1
    ## phys index mappng (julia -> actual), 1 -> -2, 2-> -1, 3 -> 0, 4 -> 1, 5 -> 2
    ## place singlets on left and bottom 
    ## tensor is projector of four virtual spin-1/2s onto a Spin=2 physical state, and \sigma_y on left and bottom legs to account for singlets
    left_idx = Index(2, "$(tag_prefix)_L")
    left_idx_dummy = Index(2, "$(tag_prefix)_Ldummy")
    bottom_idx = Index(2, "$(tag_prefix)_B") 
    bottom_idx_dummy = Index(2, "$(tag_prefix)_Bdummy")
    right_idx = Index(2, "$(tag_prefix)_R")
    top_idx = Index(2, "$(tag_prefix)_T")
    phy_idx = Index(5, "$(tag_prefix)_phys")
    akltens = ITensor(Float64, phy_idx, left_idx_dummy, bottom_idx_dummy, right_idx, top_idx)
    for iLeft in 1:2, iBottom in 1:2, iRight in 1:2, iTop in 1:2
        sum_ = iLeft + iBottom + iRight + iTop 
        if sum_ == 4 # 1111 -> 5
            akltens[5, iLeft, iBottom, iRight, iTop] = 1 
        elseif sum_ == 5 # perm(2111) -> 4
            akltens[4, iLeft, iBottom, iRight, iTop] = 1 / sqrt(4)
        elseif sum_ == 6 # perm(2211) ->  3 
            akltens[3, iLeft, iBottom, iRight, iTop] = 1 / sqrt(6)
        elseif sum_ == 7 # perm(2221) -> 2
            akltens[2, iLeft, iBottom, iRight, iTop] = 1 / sqrt(4)
        elseif sum_ == 8 # 2222 -> 1 
            akltens[1, iLeft, iBottom, iRight, iTop] = 1 
        end 
    end
    akltens = akltens * iσy(left_idx, left_idx_dummy) * iσy(bottom_idx, bottom_idx_dummy)
    return akltens
end 

function deformed_aklt_tensor(a1, a2; tag_prefix="")
    T_w = aklt_square_tensor(;tag_prefix = tag_prefix)
    idx_phys = inds(T_w)[1]
    idx_phys_new = Index(5, "$(tag_prefix)phys_deformed")
    for iLeft in 1:2, iBottom in 1:2, iRight in 1:2, iTop in 1:2
        T_w[1, iLeft, iBottom, iRight, iTop] *= a2 / sqrt(6) 
        T_w[2, iLeft, iBottom, iRight, iTop] *= a1 / sqrt(3/2)
        T_w[3, iLeft, iBottom, iRight, iTop] *= 1
        T_w[4, iLeft, iBottom, iRight, iTop] *= a1 / sqrt(3/2)
        T_w[5, iLeft, iBottom, iRight, iTop] *= a2 / sqrt(6) 
    end 
    # Dt = ITensor(idx_phys, idx_phys_new)
    # Dt[1, 1] = a2 / sqrt(6)
    # Dt[2, 2] = a1 / sqrt(3/2)
    # Dt[3, 3] = 1 
    # Dt[4, 4] = a1 / sqrt(3/2)
    # Dt[5, 5] = a2 / sqrt(6)
    return T_w # * Dt
end  

function aklt_norm_tensor(indices; a1 = sqrt(3/2), a2 = sqrt(6))
    # indices: a length-4 array of Index objects, order [L, B, R, T]
    T_Wk = deformed_aklt_tensor(a1, a2; tag_prefix = "ket")
    T_Wb = deformed_aklt_tensor(a1, a2; tag_prefix = "bra")

    # Contract the corresponding legs of ket and bra tensors
    # Result is a tensor with four legs, corresponding to the input indices
    # The contraction is over the physical indices of T_Wk and T_Wb
    function decode(idx)
        if idx == 1
            return 1, 1
        elseif idx == 2
            return 1, 2
        elseif idx == 3
            return 2, 1
        else
            return 2, 2
        end
    end

    # Create a tensor with the desired output indices
    out_tensor = ITensor(indices...)

    # For each combination of output indices, sum over the physical indices
    for iL in 1:4, iB in 1:4, iR in 1:4, iT in 1:4
        # Map each 4-d index to two 2-d indices for ket and bra
                # Map output index to ket/bra indices
                # 1->00, 2->01, 3->10, 4->11
        kL_, bL_ = decode(iL)
        kB_, bB_ = decode(iB)
        kR_, bR_ = decode(iR)
        kT_, bT_ = decode(iT)
        if (kL_ == bL_) && (kB_ == bB_)
            out_tensor[iL, iB, iR, iT] = sum([T_Wk[iPhys, kL_, kB_, kR_, kT_] * T_Wb[iPhys, bL_, bB_, bR_, bT_] for iPhys in 1:5])
        end 
    end
    return out_tensor
end 

function periodic_square_lattice(L::Int)
    g = Graph(L^2)

    # Correct row-major indexing
    idx(i, j) = (i - 1) * L + j
    
    for i in 1:L
        for j in 1:L
            current = idx(i, j)
            right = j < L ? idx(i, j + 1) : idx(i, 1)  # Wrap rightmost to first column
            down = i < L ? idx(i + 1, j) : idx(1, j)    # Wrap bottom to top row
            add_edge!(g, current, right)
            add_edge!(g, current, down)
        end
    end
    return g
end

function aklt_norm_network(L;a1 = sqrt(3/2), a2 = sqrt(6))
    inv_idx(k) = ((k - 1) ÷ L + 1, (k - 1) % L + 1)
    idx(i, j) = (i - 1) * L + j
    wrap(i) = mod(i - 1, L) + 1
    g = periodic_square_lattice(L)
    N = L^2
    indmat = [Index(4, "i$(i)j$(j)") for i in 1:N, j in 1:N]
    TN = []
    # peps = Matrix{ITensor}(undef, L, L)
    for index = 1:N 
        ii, jj = inv_idx(index)
        nbr_right  = idx(wrap(ii + 1), jj)
        nbr_left   = idx(wrap(ii - 1), jj)
        nbr_top    = idx(ii, wrap(jj + 1))
        nbr_bottom = idx(ii, wrap(jj - 1))
        indices = [indmat[min(index, nbr_left), max(index, nbr_left)], 
                   indmat[min(index, nbr_bottom), max(index, nbr_bottom)],
                   indmat[min(index, nbr_right), max(index, nbr_right)],
                   indmat[min(index, nbr_top), max(index, nbr_top)]]
        sitetensor = aklt_norm_tensor(indices; a1 = a1, a2 = a2)
        # peps[ii, jj] = sitetensor
        push!(TN, sitetensor)
    end 
    return TN
end

function dumbcontract(TN)
    s = 1 
    for t in TN
        s *= t
    end 
    return scalar(s)
end 


function ctmrg_exact_FE_density(a1,a2; χmax = 20, cutoff = 1e-8, nsteps = 400)
    sₕ = Index(4, "Right")
    sᵥ = Index(4, "Top")

    T = aklt_norm_tensor([sₕ', sᵥ', sₕ, sᵥ]; a1 = a1, a2 = a2)

    χ0 = 1
    l = Index(χ0, "Link")
    lₕ = addtags(l, "horiz")
    lᵥ = addtags(l, "vert")

    # Initial CTM
    Cₗᵤ = ITensor(lᵥ, lₕ)
    Cₗᵤ[1, 1] = 1.0
    # Initial HRTM
    Aₗ = ITensor(lᵥ, lᵥ', sₕ)
    Aₗ[lᵥ => 1, lᵥ' => 1, sₕ => 1] = 1.0

    Cₗᵤ, Aₗ = ctmrg(T, Cₗᵤ, Aₗ; χmax=χmax, cutoff=cutoff, nsteps=nsteps)

    lᵥ = commonind(Cₗᵤ, Aₗ)
    lₕ = uniqueind(Cₗᵤ, Aₗ)
    Aᵤ = replaceinds(Aₗ, lᵥ => lₕ, lᵥ' => lₕ', sₕ => sᵥ)
    ACₗ = Aₗ * Cₗᵤ * dag(Cₗᵤ')
    ACTₗ = prime(ACₗ * dag(Aᵤ') * T * Aᵤ, -1)
    exact_PF_per_site = (ACTₗ * dag(ACₗ))[]
    return log(exact_PF_per_site)
end 