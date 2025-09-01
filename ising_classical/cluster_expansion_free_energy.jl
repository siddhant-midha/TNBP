#!/usr/bin/env julia

"""
cluster_expansion_free_energy.jl

Implementation of cluster expansion for 2D Ising model free energy density using saved clusters.
Follows exactly the workflow from functions/cluster_expansion_2d_ising.jl but uses pre-computed clusters.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("../functions/Ising2D.jl")
include("../functions/BP.jl")
# include("../test/test_utils.jl") ## removing this because test utils contains duplication of the same functions used here
using Serialization
using Plots
using ITensors, Graphs

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

function load_latest_single_site_cluster_file(; size_filter="L11", weight_filter="", boundary_filter="periodic", save_dir = "../saved_clusters")
    """
    Load the most recent single-site cluster enumeration file matching criteria.
    Returns the data and filename for reference.
    """
    
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

function cluster_expansion_2d_ising_with_single_site_clusters(L::Int, β::Float64, h::Float64=0.0; 
                                                          max_weights=[4, 6, 8, 10])
    """
    Complete cluster expansion workflow for 2D Ising model using single-site saved clusters.
    Automatically loads the latest L×L single-site cluster file and handles normalization correctly.
    """
    
    println("="^80)
    println("🔥 Cluster Expansion for 2D Ising Model (Using Single-Site Clusters)")
    println("="^80)
    println("Parameters: L=$L, β=$β, h=$h")
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
    
    # Step 2: Create Ising tensor network
    println("\nStep 2: Creating Ising tensor network...")
    T = Ising2D.get_ising_tn(L, β; h=h)
    N = L^2
    println("✅ Created $N tensors for $(L)×$(L) lattice")
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    println("✅ Found $(length(edges)) edges")
    
    # Step 3: Compute BP fixed point
    println("\nStep 3: Computing BP fixed point...")
    messages = BP.get_messages(T, edges, links; random_part=0.01)
    messages = BP.message_passing(T, messages, edges, adj_mat;  α = 1., noise=0., max_iters=1000)
    println("✅ BP converged")
    
    # Get BP partition function (before normalization)
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Z_bp_full)
    println("✅ BP log partition function: $log_Z_bp_full")
    
    # Step 4: Normalize tensors
    println("\nStep 4: Normalizing tensors...")
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)
    
    # Normalization factor
    normalization_factor = sum(log.(real.(Z_l)))
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
    results["beta"] = β
    results["h"] = h
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
        # - BP contribution needs to be divided by N to get free energy density: f_BP = -log(Z_BP)/(2*N)
        # - Single-site cluster correction is already per-site, so it's the free energy density correction directly
        f_bp = -real(log_Z_bp_full) / (2 * N)  # BP free energy density
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
    
    # Step 6: Compare with exact Onsager solution
    if h == 0.0  # Onsager solution only valid for h=0
        println("\n" * "="^60)
        println("📐 Comparing with exact Onsager solution")
        println("="^60)
        
        exact_free_energy = Ising2D.free_energy(β)
        results["exact_onsager"] = exact_free_energy
        
        println("🎯 Exact Onsager free energy density: $exact_free_energy")
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
    else
        println("⚠️  Exact comparison only available for h=0 (Onsager solution)")
        results["exact_onsager"] = nothing
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
        
        if abs(phi_W) < 1e-15
            continue  # Skip negligible contributions
        end
        
        # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
        Z_W = 1.0 + 0im  # Complex number for cluster contribution
        computation_successful = true
        
        for loop_id in cluster.loop_ids
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            
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
        # The contribution to log(Z) is φ(W) * Z_W, so contribution to f is -φ(W) * Z_W / 2.0
        contribution = -phi_W * Z_W ####   / 2.0  # Factor of 2 comes from Ising model convention
        
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

function compute_single_site_cluster_correction(T_normalized, messages, edges, links, adj_mat, 
                                              cluster_data, L::Int, max_weight::Int)
    """
    Compute the cluster expansion correction using single-site cluster data.
    
    For single-site clusters: the clusters are already representatives from one site,
    and the contribution is already per-site. No need for translation multiplication.
    """
    
    # Get clusters from single-site data
    representative_clusters = cluster_data.clusters
    all_loops = cluster_data.all_loops
    
    println("📊 Using $(length(representative_clusters)) single-site clusters")
    println("   Available loops: $(length(all_loops))")
    println("   Computing per-site free energy density correction")
    
    # Filter clusters by weight
    relevant_clusters = [c for c in representative_clusters if c.weight <= max_weight]
    println("   Using $(length(relevant_clusters)) clusters (weight ≤ $max_weight)")
    
    if isempty(relevant_clusters)
        println("⚠️  No relevant clusters found, returning zero correction")
        return 0.0
    end
    
    # Compute total correction from single-site clusters
    total_correction = 0.0
    
    println("💫 Computing single-site cluster contributions...")
    
    for (i, cluster) in enumerate(relevant_clusters)
        if i % 100 == 0 || i <= 10
            println("  Processing cluster $i/$(length(relevant_clusters))...")
        end
        
        # Compute Ursell function φ(W)
        phi_W = ursell_function(cluster)
        
        if abs(phi_W) < 1e-15
            continue  # Skip negligible contributions
        end
        
        # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
        Z_W = 1.0 + 0im  # Complex number for cluster contribution
        
        for loop_id in cluster.loop_ids
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            
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
                continue
            end
        end
        
        # Add weighted contribution: φ(W) * Z_W
        # For free energy: f = -log(Z), so correction is φ(W) * Z_W, but we want the contribution to f
        # The contribution to log(Z) is φ(W) * Z_W, so contribution to f is -φ(W) * Z_W / (2*N)
        # But since we're computing per-site correction, we don't divide by N
        contribution = -phi_W * Z_W / 2.0  # Factor of 2 comes from Ising model convention
        total_correction += contribution
        
        if abs(contribution) > 1e-10 && i <= 20
            println("    Cluster $i: φ=$phi_W, Z_W=$Z_W, f contribution=$contribution")
        end
    end
    
    println("✅ Single-site cluster correction (free energy density): $total_correction")
    return real(total_correction)  # Return free energy density correction per site
end

function cluster_expansion_2d_ising_with_saved_clusters(cluster_file::String, L::Int, β::Float64, h::Float64=0.0; 
                                                       max_weights=[4, 6, 8, 10])
    """
    Complete cluster expansion workflow for 2D Ising model using saved clusters.
    Follows exactly the pattern from functions/cluster_expansion_2d_ising.jl
    """
    
    println("="^80)
    println("🔥 Cluster Expansion for 2D Ising Model (Using Saved Clusters)")
    println("="^80)
    println("Parameters: L=$L, β=$β, h=$h")
    println("Cluster weights: $max_weights")
    println()
    
    # Step 1: Create Ising tensor network
    println("Step 1: Creating Ising tensor network...")
    T = Ising2D.get_ising_tn(L, β; h=h)
    N = L^2
    println("✅ Created $N tensors for $(L)×$(L) lattice")
    
    # Get adjacency structure
    adj_mat, edges, links = BP.get_adj_mat(T)
    println("✅ Found $(length(edges)) edges")
    
    # Step 2: Compute BP fixed point
    println("\nStep 2: Computing BP fixed point...")
    messages = BP.get_messages(T, edges, links;random_part=0.01)
    messages = BP.message_passing(T, messages, edges, adj_mat; max_iters=1000)
    println("✅ BP converged")
    
    # Get BP partition function (before normalization)
    Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
    log_Z_bp_full = log(Z_bp_full)
    println("✅ BP log partition function: $log_Z_bp_full")
    
    # Step 3: Normalize tensors
    println("\nStep 3: Normalizing tensors...")
    Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
    T_normalized = BP.normalize_tensors(T, Z_l)
    
    # Normalization factor
    normalization_factor = sum(log.(real.(Z_l)))
    println("✅ Normalization factor: $normalization_factor")
    
    # Verify normalization (BP on normalized tensors should give 1)
    Z_bp_normalized = BP.mean_free_partition_fn(1:N, T_normalized, messages, adj_mat)
    println("✅ Normalized BP partition function: $Z_bp_normalized (should ≈ 1)")
    
    # Step 4: Load saved clusters and compute corrections
    println("\nStep 4: Loading saved clusters...")
    cluster_data, cluster_summary = load_cluster_data(cluster_file)
    
    results = Dict{String, Any}()
    results["L"] = L
    results["beta"] = β
    results["h"] = h
    results["bp_log_Z"] = log_Z_bp_full  # Keep as raw log partition function
    results["normalization_factor"] = normalization_factor
    results["cluster_corrections"] = Dict()
    
    for max_weight in max_weights
        println("\n" * "="^60)
        println("🔍 Computing cluster expansion up to weight $max_weight")
        println("="^60)
        
        correction = compute_cluster_correction_with_saved_data(
            T_normalized, messages, edges, links, adj_mat, 
            cluster_data, L, max_weight)
        
        # Total log partition function with cluster correction
        # According to polymer theory: log(Z) = log(Z_BP) + Σ φ(W) Z_W
        # For free energy density: f = -log(Z)/(2*N), so correction is added directly
        log_Z_corrected = log_Z_bp_full + correction * N  # correction is per-site, multiply by N for total log(Z)
        
        results["cluster_corrections"][max_weight] = Dict(
            "correction" => correction,  # This is the free energy density correction
            "total_log_Z" => log_Z_corrected,
            "improvement" => correction
        )
        
        println("📊 Results for weight $max_weight:")
        println("  Cluster correction: $correction")
        println("  Total log partition function: $log_Z_corrected")
    end
    
    # Step 5: Compare with exact Onsager solution
    if h == 0.0  # Onsager solution only valid for h=0
        println("\n" * "="^60)
        println("📐 Comparing with exact Onsager solution")
        println("="^60)
        
        exact_free_energy = Ising2D.free_energy(β)
        exact_log_Z = -exact_free_energy * (2 * N)  # Convert free energy back to log Z
        results["exact_onsager"] = exact_free_energy
        results["exact_log_Z"] = exact_log_Z
        
        println("🎯 Exact Onsager free energy: $exact_free_energy")
        println("🎯 Exact log partition function: $exact_log_Z")
        println("\n📈 Comparison with exact result:")
        println("  BP approximation error: $(abs(results["bp_log_Z"] - exact_log_Z))")
        
        for max_weight in max_weights
            corrected_log_Z = results["cluster_corrections"][max_weight]["total_log_Z"]
            error = abs(corrected_log_Z - exact_log_Z)
            improvement = abs(results["bp_log_Z"] - exact_log_Z) - error
            
            println("  Weight $max_weight error: $error (improvement: $improvement)")
        end
    else
        println("⚠️  Exact comparison only available for h=0 (Onsager solution)")
        results["exact_onsager"] = nothing
        results["exact_log_Z"] = nothing
    end
    
    return results
end

function compute_cluster_correction_with_saved_data(T_normalized, messages, edges, links, adj_mat, 
                                                   cluster_data, L::Int, max_weight::Int)
    """
    Compute the cluster expansion correction using saved cluster data.
    
    CRITICAL: The saved clusters are "representatives" from one site with translational 
    symmetry removed. For periodic boundary conditions, we need to multiply by the 
    number of equivalent translations (N = L²) to get the total contribution.
    """
    
    # Get clusters from saved data
    if hasfield(typeof(cluster_data), :clusters_by_site)
        # All-sites enumeration data - use site 1 clusters (PBC makes all equivalent)
        representative_clusters = cluster_data.clusters_by_site[1]
    else
        # Single-site enumeration data
        representative_clusters = cluster_data.clusters
    end
    
    all_loops = cluster_data.all_loops
    N = L^2  # Total number of sites
    
    println("📊 Using $(length(representative_clusters)) representative clusters from saved data")
    println("   Available loops: $(length(all_loops))")
    println("   Will multiply by N=$N for translational equivalence in PBC")
    
    # Filter clusters by weight
    relevant_clusters = [c for c in representative_clusters if c.weight <= max_weight]
    println("   Using $(length(relevant_clusters)) representative clusters (weight ≤ $max_weight)")
    
    if isempty(relevant_clusters)
        println("⚠️  No relevant clusters found, returning zero correction")
        return 0.0
    end
    
    # Compute total correction from representatives
    representative_correction = 0.0
    
    println("💫 Computing cluster contributions...")
    # Process all relevant clusters (no arbitrary limit like in demo)
    
    for (i, cluster) in enumerate(relevant_clusters)
        if i % 100 == 0 || i <= 10
            println("  Processing representative cluster $i/$(length(relevant_clusters))...")
        end
        
        # Compute Ursell function φ(W)
        phi_W = ursell_function(cluster)
        
        if abs(phi_W) < 1e-15
            continue  # Skip negligible contributions
        end
        
        # Compute cluster correction Z_W = ∏_i Z_{l_i}^{η_i}
        Z_W = 1.0 + 0im  # Complex number for cluster contribution
        
        for loop_id in cluster.loop_ids
            multiplicity = cluster.multiplicities[loop_id]
            loop = all_loops[loop_id]
            
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
                continue
            end
        end
        
        # Add weighted contribution: φ(W) * Z_W
        contribution = phi_W * Z_W
        representative_correction += contribution
        
        if abs(contribution) > 1e-10 && i <= 20
            println("    Representative cluster $i: φ=$phi_W, Z_W=$Z_W, contribution=$contribution")
        end
    end
    
    # CRITICAL: For PBC, multiply by N to account for all translational copies
    # But since we're computing FREE ENERGY DENSITY, this factor of N cancels out
    # when we divide by N in the final free energy density calculation.
    # So the total correction to the LOG partition function is N * representative_correction,
    # but the correction to the FREE ENERGY DENSITY is just representative_correction.
    total_log_Z_correction = N * representative_correction
    free_energy_density_correction = representative_correction  # This is what gets returned
    
    println("✅ Representative cluster correction: $representative_correction")
    println("✅ Total log(Z) correction: $total_log_Z_correction")
    println("✅ Free energy density correction: $free_energy_density_correction")
    return real(free_energy_density_correction)  # Return free energy density correction
end

# Plot free energy vs beta
function plot_free_energy_vs_beta_with_saved_clusters(cluster_file::String, β_range::Vector{Float64}, 
                                                     max_weights::Vector{Int}=[4,6,8,10]; save_path::String="")
    """
    Plot free energy density as a function of β using saved clusters.
    """
    
    println("📊 Creating free energy vs β plot using saved clusters")
    println("   β range: $(β_range[1]) to $(β_range[end]) ($(length(β_range)) points)")
    println("   Weight truncations: $max_weights")
    
    # Load cluster data to get lattice size
    cluster_data, cluster_summary = load_cluster_data(cluster_file)
    L = cluster_summary["lattice_size"]
    
    # Initialize plot
    p = plot(title="2D Ising Free Energy Density vs β (Cluster Expansion)", 
             xlabel="β", ylabel="Free Energy Density f",
             legend=:bottomleft, size=(800, 600),
             dpi=300)
    
    # Compute exact Onsager solution
    println("\n🎯 Computing exact Onsager solution...")
    exact_f = [Ising2D.free_energy(β) for β in β_range]
    plot!(p, β_range, exact_f, linewidth=3, color=:black, 
          label="Exact (Onsager)", linestyle=:solid)
    
    # Compute BP baseline (more efficiently - just compute BP without clusters)
    println("\n📊 Computing BP baseline...")
    bp_f = Float64[]
    for (j, β) in enumerate(β_range)
        if j % 5 == 1
            println("   BP β = $β ($(j)/$(length(β_range)))")
        end
        
        try
            # Quick BP computation without clusters
            T = Ising2D.get_ising_tn(L, β; h=0.0)
            N = L^2
            adj_mat, edges, links = BP.get_adj_mat(T)
            messages = BP.get_messages(T, edges, links;random_part=0.01)
            messages = BP.message_passing(T, messages, edges, adj_mat; max_iters=1000)
            Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
            
            # BP free energy density: f = -log(Z_BP) / (2*N)
            bp_log_Z = log(Z_bp_full)
            f = -real(bp_log_Z) / (2 * L^2)  # Divide by N for BP
            push!(bp_f, f)
        catch e
            println("     ⚠️ Error at β=$β: $e")
            push!(bp_f, NaN)
        end
    end
    
    # Plot BP baseline
    plot!(p, β_range, bp_f, linewidth=2, color=:gray, 
          label="BP (no corrections)", linestyle=:dot)
    
    # Colors for different weights
    colors = [:red, :blue, :green, :orange, :purple, :brown]
    
    # Compute cluster expansion for each weight
    for (i, max_weight) in enumerate(max_weights)
        println("\n📈 Computing cluster expansion for weight $max_weight...")
        
        cluster_f = Float64[]
        
        for (j, β) in enumerate(β_range)
            if j % 5 == 1
                println("   β = $β ($(j)/$(length(β_range)))")
            end
            
            try
                results = cluster_expansion_2d_ising_with_saved_clusters(cluster_file, L, β, 0.0; max_weights=[max_weight])
                # For cluster expansion: f = -log(Z_BP) / (2*N) + correction
                # where correction is already the free energy density correction
                bp_log_Z = results["bp_log_Z"]
                correction = results["cluster_corrections"][max_weight]["correction"]
                # Both BP and correction are now free energy density terms
                f = -real(bp_log_Z) / (2 * L^2) + real(correction)
                push!(cluster_f, f)
            catch e
                println("     ⚠️ Error at β=$β: $e")
                push!(cluster_f, NaN)
            end
        end
        
        # Plot this weight truncation
        color = colors[min(i, length(colors))]
        plot!(p, β_range, cluster_f, linewidth=2, color=color,
              label="Weight ≤ $max_weight", linestyle=:dash)
    end
    
    # Save plot
    if !isempty(save_path)
        println("\n💾 Saving plot to: $save_path")
        savefig(p, save_path)
    end
    
    return p
end

# Plot free energy vs beta using single-site clusters (EFFICIENT VERSION)
function plot_free_energy_vs_beta_with_single_site_clusters_efficient(L::Int, β_range::Vector{Float64}, 
                                                                      max_weights::Vector{Int}=[4,6,8,10]; save_path::String="")
    """
    Plot free energy density as a function of β using single-site clusters.
    EFFICIENT: Loads cluster data once and computes contributions once per β value.
    """
    
    println("📊 Creating free energy vs β plot using single-site clusters (EFFICIENT)")
    println("   Lattice size: $(L)×$(L)")
    println("   β range: $(β_range[1]) to $(β_range[end]) ($(length(β_range)) points)")
    println("   Weight truncations: $max_weights")
    
    # Load single-site cluster data once
    println("\n📂 Loading single-site cluster data once...")
    cluster_data = nothing
    cluster_filename = ""
    
    try
        # Try to load any file that has weight >= maximum(max_weights)
        min_required_weight = maximum(max_weights)
        println("⚠️  Looking for clusters with weight >= $min_required_weight")
        
        # First try to load any weight that meets our requirement
        try
            cluster_data, cluster_filename = load_latest_single_site_cluster_file(size_filter="L$L", weight_filter="w$(min_required_weight)")
            
            # Check if the loaded file has sufficient weight
            if cluster_data.max_weight >= min_required_weight
                println("✅ Loaded single-site cluster data from: $(cluster_filename)")
                println("   File max weight: $(cluster_data.max_weight) (>= required: $min_required_weight)")
            else
                println("⚠️  Loaded file has max weight $(cluster_data.max_weight) < required $min_required_weight")
                println("   Will still proceed but some weight truncations may be incomplete")
            end
        catch e
            error("❌ Could not load any single-site cluster file for L=$L: $e")
        end
    catch e
        error("❌ Could not load single-site cluster data for L=$L: $e")
    end
    
    # Initialize plot
    p = plot(title="2D Ising Free Energy Density vs β (Single-Site Cluster Expansion)", 
             xlabel="β", ylabel="Free Energy Density f",
             legend=:bottomleft, size=(800, 600),
             dpi=300)
    
    # Compute exact Onsager solution
    println("\n🎯 Computing exact Onsager solution...")
    exact_f = [Ising2D.free_energy(β) for β in β_range]
    plot!(p, β_range, exact_f, linewidth=3, color=:black, 
          label="Exact (Onsager)", linestyle=:solid)
    
    # Initialize storage for all β points
    bp_f = Float64[]
    cluster_f_dict = Dict{Int, Vector{Float64}}()
    for max_weight in max_weights
        cluster_f_dict[max_weight] = Float64[]
    end
    
    max_weight_overall = maximum(max_weights)
    
    # Compute for each β value
    println("\n📊 Computing for each β value...")
    for (j, β) in enumerate(β_range)
        if j % 5 == 1 || length(β_range) <= 10
            println("   Processing β = $β ($(j)/$(length(β_range)))")
        end
        
        try
            # Create Ising tensor network
            T = Ising2D.get_ising_tn(L, β; h=0.0)
            N = L^2
            adj_mat, edges, links = BP.get_adj_mat(T)
            messages = BP.get_messages(T, edges, links;random_part=0.03)
            messages = BP.message_passing(T, messages, edges, adj_mat; max_iters=1000)
            
            Z_l = BP.get_fixed_point_list(T, messages, adj_mat)
            log_Z_sum = sum(log.(real.(Z_l)))

            # Get BP free energy density
            Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
            bp_log_Z = log(Z_bp_full)
            f_bp = - log_Z_sum / (2 * L^2)
            push!(bp_f, f_bp)
            
            # Normalize tensors
            
            T_normalized = BP.normalize_tensors(T, Z_l)
            
            # Compute ALL cluster contributions once for this β
            cluster_contributions = compute_all_single_site_cluster_contributions(
                T_normalized, messages, edges, links, adj_mat, 
                cluster_data, L, max_weight_overall)
            
            # Sum contributions for different weight truncations
            for max_weight in max_weights
                correction = sum(contrib["contribution"] for contrib in cluster_contributions if contrib["weight"] <= max_weight)
                f_corrected = f_bp + correction
                push!(cluster_f_dict[max_weight], f_corrected)
            end
            
        catch e
            println("     ⚠️ Error at β=$β: $e")
            push!(bp_f, NaN)
            for max_weight in max_weights
                push!(cluster_f_dict[max_weight], NaN)
            end
        end
    end
    
    # Plot BP baseline
    plot!(p, β_range, bp_f, linewidth=2, color=:gray, 
          label="BP (no corrections)", linestyle=:dot)
    
    # Colors for different weights
    colors = [:red, :blue, :green, :orange, :purple, :brown]
    
    # Plot cluster expansion for each weight
    for (i, max_weight) in enumerate(max_weights)
        color = colors[min(i, length(colors))]
        plot!(p, β_range, cluster_f_dict[max_weight], linewidth=2, color=color,
              label="Weight ≤ $max_weight", linestyle=:dash)
    end
    
    # Save plot
    if !isempty(save_path)
        println("\n💾 Saving plot to: $save_path")
        savefig(p, save_path)
    end
    
    # Create error plot showing absolute differences from exact solution
    println("\n📊 Creating error plot showing absolute differences from exact solution...")
    p_error = plot(title="Absolute Error vs β (Distance from Exact Onsager Solution)", 
                   xlabel="β", ylabel="Absolute Error |f - f_exact|",
                   legend=:topright, size=(800, 600),
                   dpi=300, yscale=:log10)
    
    # Plot BP error
    bp_error = abs.(bp_f .- exact_f)
    plot!(p_error, β_range, bp_error, linewidth=2, color=:gray, 
          label="BP (no corrections)", linestyle=:dot)
    
    # Plot cluster expansion errors for each weight
    for (i, max_weight) in enumerate(max_weights)
        cluster_error = abs.(cluster_f_dict[max_weight] .- exact_f)
        color = colors[min(i, length(colors))]
        plot!(p_error, β_range, cluster_error, linewidth=2, color=color,
              label="Weight ≤ $max_weight", linestyle=:dash)
    end
    
    # Save error plot
    if !isempty(save_path)
        error_save_path = replace(save_path, ".png" => "_errors.png")
        println("💾 Saving error plot to: $error_save_path")
        savefig(p_error, error_save_path)
    end
    
    return p, p_error
end

# Plot free energy vs beta using single-site clusters
function plot_free_energy_vs_beta_with_single_site_clusters(L::Int, β_range::Vector{Float64}, 
                                                          max_weights::Vector{Int}=[4,6,8,10]; save_path::String="")
    """
    Plot free energy density as a function of β using single-site clusters.
    Automatically loads the latest single-site cluster data.
    """
    
    println("📊 Creating free energy vs β plot using single-site clusters")
    println("   Lattice size: $(L)×$(L)")
    println("   β range: $(β_range[1]) to $(β_range[end]) ($(length(β_range)) points)")
    println("   Weight truncations: $max_weights")
    
    # Initialize plot
    p = plot(title="2D Ising Free Energy Density vs β (Single-Site Cluster Expansion)", 
             xlabel="β", ylabel="Free Energy Density f",
             legend=:bottomleft, size=(800, 600),
             dpi=300)
    
    # Compute exact Onsager solution
    println("\n🎯 Computing exact Onsager solution...")
    exact_f = [Ising2D.free_energy(β) for β in β_range]
    plot!(p, β_range, exact_f, linewidth=3, color=:black, 
          label="Exact (Onsager)", linestyle=:solid)
    
    # Compute BP baseline (more efficiently - just compute BP without clusters)
    println("\n📊 Computing BP baseline...")
    bp_f = Float64[]
    for (j, β) in enumerate(β_range)
        if j % 5 == 1
            println("   BP β = $β ($(j)/$(length(β_range)))")
        end
        
        try
            # Quick BP computation without clusters
            T = Ising2D.get_ising_tn(L, β; h=0.0)
            N = L^2
            adj_mat, edges, links = BP.get_adj_mat(T)
            messages = BP.get_messages(T, edges, links;random_part=0.01)
            messages = BP.message_passing(T, messages, edges, adj_mat; max_iters=1000)
            Z_bp_full = BP.mean_free_partition_fn(1:N, T, messages, adj_mat)
            
            # BP free energy density: f = -log(Z_BP) / (2*N)
            bp_log_Z = log(Z_bp_full)
            f = -real(bp_log_Z) / (2 * L^2)  # Divide by N for BP
            push!(bp_f, f)
        catch e
            println("     ⚠️ Error at β=$β: $e")
            push!(bp_f, NaN)
        end
    end
    
    # Plot BP baseline
    plot!(p, β_range, bp_f, linewidth=2, color=:gray, 
          label="BP (no corrections)", linestyle=:dot)
    
    # Colors for different weights
    colors = [:red, :blue, :green, :orange, :purple, :brown]
    
    # Compute cluster expansion for each weight
    for (i, max_weight) in enumerate(max_weights)
        println("\n📈 Computing single-site cluster expansion for weight $max_weight...")
        
        cluster_f = Float64[]
        
        for (j, β) in enumerate(β_range)
            if j % 5 == 1
                println("   β = $β ($(j)/$(length(β_range)))")
            end
            
            try
                results = cluster_expansion_2d_ising_with_single_site_clusters(L, β, 0.0; max_weights=[max_weight])
                # Get corrected free energy density directly
                f_corrected = results["cluster_corrections"][max_weight]["f_corrected"]
                push!(cluster_f, f_corrected)
            catch e
                println("     ⚠️ Error at β=$β: $e")
                push!(cluster_f, NaN)
            end
        end
        
        # Plot this weight truncation
        color = colors[min(i, length(colors))]
        plot!(p, β_range, cluster_f, linewidth=2, color=color,
              label="Weight ≤ $max_weight", linestyle=:dash)
    end
    
    # Save plot
    if !isempty(save_path)
        println("\n💾 Saving plot to: $save_path")
        savefig(p, save_path)
    end
    
    return p
end

# Main execution function
function main()
    println("🚀 Cluster Expansion Free Energy Calculator (Using Single-Site Clusters - EFFICIENT)")
    println("="^80)
    
    # Parameters - automatically use L10 since we have single-site data for it
    L = 6
    β_range = collect(0.36:0.01:0.38)  # β from 0.1 to 1.6 with fewer points for faster computation
    max_weights = [10]  # Use available weights for L10 (we have w8 and w9)
    
    println("📋 Configuration:")
    println("   Lattice size: $(L)×$(L)")
    println("   β range: $(β_range[1]) to $(β_range[end]) ($(length(β_range)) points)")
    println("   Weight truncations: $max_weights")
    println("   Using EFFICIENT computation: load once, compute once per β")
    
    # Create plot using efficient single-site clusters
    save_path = "../visualization/cluster_expansion_free_energy_L$(L)_single_site_efficient.png"
    plot_obj, error_plot_obj = plot_free_energy_vs_beta_with_single_site_clusters_efficient(L, β_range, max_weights; save_path=save_path)
    
    # Display final results
    println("\n✅ Analysis complete!")
    println("📊 Free energy plot saved to: $save_path")
    error_save_path = replace(save_path, ".png" => "_errors.png")
    println("📊 Error plot saved to: $error_save_path")
    
    return plot_obj, error_plot_obj
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end