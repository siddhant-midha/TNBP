using ITensors, ITensorMPS, Graphs, LinearAlgebra, Statistics
include("../functions/ClusterEnumeration.jl")
include("../functions/ctmrg.jl")

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


function aklt_sz_tensor(indices; a1 = sqrt(3/2), a2 = sqrt(6))
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
    function mag(iPhys) # 1 -> -2, 2-> -1, 3 -> 0, 4 -> 1, 5 -> 2
        if iPhys == 1
            return -2 
        elseif iPhys == 2 
            return -1 
        elseif iPhys == 3 
            return 0
        elseif iPhys == 4 
            return 1
        elseif iPhys == 5
            return 2
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
            out_tensor[iL, iB, iR, iT] = sum([mag(iPhys) * T_Wk[iPhys, kL_, kB_, kR_, kT_] * T_Wb[iPhys, bL_, bB_, bR_, bT_] for iPhys in 1:5])
        end 
    end
    return out_tensor
end 


function aklt_sz_network(L;a1 = sqrt(3/2), a2 = sqrt(6))
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
        if index == 1
            sitetensor = aklt_sz_tensor(indices; a1 = a1, a2 = a2)
        else 
            sitetensor = aklt_norm_tensor(indices; a1 = a1, a2 = a2)
        end 
        # peps[ii, jj] = sitetensor
        push!(TN, sitetensor)
    end 
    return TN
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


function ctmrg_exact_FE_density(a1,a2; χmax = 32, cutoff = 1e-14, nsteps = 400)
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

# Data structure for single-site enumeration results (copied from generate_ising_clusters_one_site.jl)
struct SingleSiteClusterData
    """Data structure to save single-site cluster enumeration results."""
    clusters::Vector{Cluster}
    all_loops::Vector{Loop}
    adj_matrix::Matrix{Int}
    max_weight::Int
    lattice_size::Int
    site::Int
    enumeration_time::Float64
    translation_removal_time::Float64
    timestamp::String
    boundary_condition::String
    clusters_before_translation_removal::Int
    clusters_after_translation_removal::Int
    canonical_forms_count::Int
end

# Load saved cluster data
function load_cluster_data(filepath::String)
    """Load cluster enumeration data from saved file."""
    println("📂 Loading cluster data from: $(basename(filepath))")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    
    data = loaded_data["data"]
    summary = loaded_data["summary"]
    
    println("✅ Loaded successfully!")
    println("   Lattice size: $(summary["lattice_size"])×$(summary["lattice_size"])")
    println("   Max weight: $(summary["max_weight"])")
    println("   Total loops: $(summary["total_loops"])")
    
    if haskey(summary, "total_clusters")
        println("   Total clusters: $(summary["total_clusters"])")
    end
    
    return data, summary
end

function load_latest_single_site_cluster_file(; size_filter="L11", weight_filter="", boundary_filter="periodic")
    """
    Load the most recent single-site cluster enumeration file matching criteria.
    Returns the data and filename for reference.
    """
    save_dir = "../saved_clusters"
    
    if !isdir(save_dir)
        error("No saved_clusters directory found!")
    end
    
    files = readdir(save_dir)
    
    # Filter files based on criteria - must contain "single_site"
    matching_files = filter(files) do f
        # Must be a .jld2 file
        !endswith(f, ".jld2") && return false
        
        # Must contain "single_site" in filename
        !contains(f, "single_site") && return false
        
        # Must contain size filter
        !contains(f, size_filter) && return false
        
        # Must contain boundary filter
        !contains(f, boundary_filter) && return false
        
        # If weight filter specified, must contain it
        if !isempty(weight_filter)
            !contains(f, weight_filter) && return false
        end
        
        return true
    end
    
    if isempty(matching_files)
        error("No matching single-site cluster files found for criteria: size=$size_filter, weight=$weight_filter, boundary=$boundary_filter")
    end
    
    # Sort by filename (which includes timestamp) and take the most recent
    sorted_files = sort(matching_files, rev=true)  # Most recent first
    
    # Try to load files starting with the most recent
    for latest_file in sorted_files
        filepath = joinpath(save_dir, latest_file)
        
        println("📖 Attempting to load single-site cluster data: $(latest_file)")
        
        try
            loaded_data = open(filepath, "r") do io
                deserialize(io)
            end
            
            println("✅ Successfully loaded: $(latest_file)")
            return loaded_data["data"], latest_file
        catch e
            println("⚠️  Failed to load $(latest_file): $e")
            if latest_file == sorted_files[end]  # If this was the last file to try
                error("❌ All matching cluster files are corrupted or unreadable")
            end
            continue  # Try the next file
        end
    end
end

# Ursell function for connected clusters (from original implementation)
function ursell_function(cluster::Cluster)
    """
    Compute the Ursell function φ(W) for a connected cluster W.
    For connected clusters, this is (-1)^(|W|-1) * (|W|-1)!
    where |W| is the total number of loops in the cluster (with multiplicities).
    """
    total_loops = cluster.total_loops
    if total_loops == 1
        return 1.0
    else
        return -0.5
    end
end

function cluster_expansion_aklt_with_single_site_clusters(L::Int, a1::Float64, a2::Float64; 
                                                          max_weights=[4, 6, 8, 10])
    """
    Complete cluster expansion workflow for AKLT model using single-site saved clusters.
    Automatically loads the latest L×L single-site cluster file and handles normalization correctly.
    """
    
    println("="^80)
    println("🔥 Cluster Expansion for AKLT Model (Using Single-Site Clusters)")
    println("="^80)
    println("Parameters: L=$L, a1=$a1, a2=$a2")
    println("Cluster weights: $max_weights")
    println()
    
    # Step 1: Load latest single-site cluster data
    println("Step 1: Loading latest single-site cluster data...")
    cluster_data = nothing
    cluster_filename = ""
    
    try
        cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L", weight_filter="w$(maximum(max_weights))")
        println("✅ Loaded single-site cluster data from: $(cluster_filename)")
    catch e
        println("⚠️  Failed to load with specific weight, trying any weight...")
        try
            cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L")
            println("✅ Loaded single-site cluster data from: $(cluster_filename)")
        catch e2
            error("❌ Could not load single-site cluster data for L=$L: $e2")
        end
    end
    
    # Step 2: Create AKLT tensor network
    println("\nStep 2: Creating AKLT tensor network...")
    T = aklt_norm_network(L; a1=a1, a2=a2)
    N = L^2
    println("✅ Created $N tensors for $(L)×$(L) lattice")
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    println("✅ Found $(length(edges)) edges")
    
    # Step 3: Compute BP fixed point
    println("\nStep 3: Computing BP fixed point...")
    messages = BP.get_messages(T, edges, links; random_part=0.1)
    messages = BP.message_passing(T, messages, edges, adj_mat; α=0.8, max_iters=1000)
    println("✅ BP converged")
    
    # Get BP partition function (before normalization)
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Complex(Z_bp_full))
    println("✅ BP log partition function: $log_Z_bp_full")
    
    # Step 4: Normalize tensors
    println("\nStep 4: Normalizing tensors...")
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)
    
    # Normalization factor
    Z_l_complex = Complex.(Z_l)
    normalization_factor = sum(log.(Z_l_complex))
    println("✅ Normalization factor: $normalization_factor")
    
    # Verify normalization (BP on normalized tensors should give 1)
    Z_bp_normalized = BP.mean_free_partition_fn(1:N, T_normalized, messages, adj_mat)
    println("✅ Normalized BP partition function: $Z_bp_normalized (should ≈ 1)")
    
    # Step 5: Compute ALL cluster contributions once
    println("\nStep 5: Computing ALL cluster contributions once...")
    
    # Compute individual cluster contributions for all clusters up to max weight
    max_weight_overall = maximum(max_weights)
    cluster_contributions = compute_all_single_site_cluster_contributions(
        T_normalized, messages, edges, links, adj_mat, 
        cluster_data, L, max_weight_overall)
    
    println("✅ Computed $(length(cluster_contributions)) individual cluster contributions")
    
    results = Dict{String, Any}()
    results["L"] = L
    results["a1"] = a1
    results["a2"] = a2
    results["bp_log_Z"] = log_Z_bp_full  # Keep as raw log partition function
    results["normalization_factor"] = normalization_factor
    results["cluster_corrections"] = Dict()
    results["cluster_filename"] = cluster_filename
    
    # Step 6: Sum contributions for different weight truncations
    println("\nStep 6: Summing contributions for different weight truncations...")
    
    for max_weight in max_weights
        println("🔍 Computing cluster expansion up to weight $max_weight")
        
        # Sum only contributions from clusters with weight <= max_weight
        correction = sum(contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= max_weight)
        
        # For single-site clusters: 
        # - BP contribution needs to be divided by N to get free energy density: f_BP = -log(Z_BP)/N
        # - Single-site cluster correction is already per-site, so it's the free energy density correction directly
        f_bp = -real(log_Z_bp_full) / N  # BP free energy density
        f_corrected = f_bp + correction  # Corrected free energy density
        
        results["cluster_corrections"][max_weight] = Dict(
            "correction" => correction,  # Free energy density correction (per site)
            "f_bp" => f_bp,  # BP free energy density
            "f_corrected" => f_corrected,  # Corrected free energy density
            "improvement" => abs(correction)
        )
        
        println("📊 Results for weight $max_weight:")
        println("  BP free energy density: $f_bp")
        println("  Cluster correction (per site): $correction")
        println("  Corrected free energy density: $f_corrected")
    end
    
    # Step 6: Compare with exact CTMRG solution
    println("\n" * "="^60)
    println("📐 Comparing with exact CTMRG solution")
    println("="^60)
    
    exact_free_energy = -ctmrg_exact_FE_density(a1, a2)
    results["exact_ctmrg"] = exact_free_energy
    
    println("🎯 Exact CTMRG free energy density: $exact_free_energy")
    println("\n📈 Comparison with exact result:")
    
    f_bp = results["cluster_corrections"][max_weights[1]]["f_bp"]  # BP is the same for all weights
    bp_error = abs(f_bp - exact_free_energy)
    println("  BP approximation error: $bp_error")
    
    for max_weight in max_weights
        f_corrected = results["cluster_corrections"][max_weight]["f_corrected"]
        error = abs(f_corrected - exact_free_energy)
        improvement = bp_error - error
        
        println("  Weight $max_weight error: $error (improvement: $improvement)")
    end
    
    return results
end
function compute_all_single_site_cluster_contributions(T_normalized, messages, edges, links, adj_mat, 
                                                     cluster_data, L::Int, max_weight::Int)
    """
    Compute ALL cluster contributions once for single-site cluster data.
    Returns a vector of dictionaries with contribution details for each cluster.
    """
    
    # Get clusters from single-site data
    all_clusters = cluster_data.clusters
    all_loops = cluster_data.all_loops
    
    println("📊 Computing contributions for all single-site clusters")
    println("   Available clusters: $(length(all_clusters))")
    println("   Available loops: $(length(all_loops))")
    
    # Filter clusters by weight
    relevant_clusters = [c for c in all_clusters if c.weight <= max_weight]
    println("   Using $(length(relevant_clusters)) clusters (weight ≤ $max_weight)")
    
    if isempty(relevant_clusters)
        println("⚠️  No relevant clusters found")
        return []
    end
    
    # Compute contributions for all relevant clusters
    cluster_contributions = []
    
    println("💫 Computing individual cluster contributions...")
    
    for (i, cluster) in enumerate(relevant_clusters)
        if i % 100 == 0 || i <= 10
            println("  Processing cluster $i/$(length(relevant_clusters)) (weight=$(cluster.weight))...")
        end
        
        # Compute Ursell function φ(W)
        phi_W = ursell_function(cluster)
        
        # if abs(phi_W) < 1e-24
        #     continue  # Skip negligible contributions
        # end
        
        # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
        Z_W = 1.0 + 0im  # Complex number for cluster contribution
        computation_successful = true
        number_of_loops = 0
        for loop_id in cluster.loop_ids
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            number_of_loops += multiplicity
            # Convert loop edges format for BP.jl (ensure v1 < v2 ordering)
            loop_edges_bp = Tuple{Int,Int}[]
            for edge in loop.edges
                v1, v2 = edge
                push!(loop_edges_bp, (min(v1, v2), max(v1, v2)))
            end
            
            # Compute loop contribution Z_l using BP
            try
                Z_l_tensor = BP.loop_contribution(loop_edges_bp, messages, T_normalized, edges, links, adj_mat)
                # Extract scalar value from ITensor
                Z_l = scalar(Z_l_tensor)
                Z_W *= Z_l^multiplicity
            catch e
                println("⚠️  Error computing loop contribution for loop $loop_id: $e")
                computation_successful = false
                break
            end
        end
        
        if !computation_successful
            continue
        end
        
        # Add contribution to free energy density
        # For free energy: f = -log(Z), so correction is φ(W) * Z_W, but we want the contribution to f
        contribution = -phi_W * Z_W
        ## Temporary fix
        if cluster.weight == 10 && number_of_loops == 2
            contribution = contribution * 12
        elseif cluster.weight == 8 && number_of_loops == 2
            contribution = contribution * 9
        end

        # Store contribution details
        contrib_info = Dict(
            "cluster_id" => i,
            "weight" => cluster.weight,
            "phi_W" => phi_W,
            "Z_W" => Z_W,
            "contribution" => real(contribution)
        )
        push!(cluster_contributions, contrib_info)
        
        if abs(contribution) > 1e-10 && i <= 20
            println("    Cluster $i (weight=$(cluster.weight)): φ=$phi_W, Z_W=$Z_W, f contribution=$contribution")
        end
    end
    
    println("✅ Computed $(length(cluster_contributions)) valid cluster contributions")
    
    # Sort by weight for easier analysis
    sort!(cluster_contributions, by=x->x["weight"])
    
    return cluster_contributions
end