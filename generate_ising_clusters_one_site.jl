#!/usr/bin/env julia

"""
generate_ising_clusters_one_site.jl

Generate and save connected cluster enumerations for a single site on square lattice Ising models.
Supports both periodic and open boundary conditions.
Includes translational symmetry removal - only keeps one representative from each equivalence class.

Usage:
    julia generate_ising_clusters_one_site.jl --size 6 --weight 4 --site 1 --boundary periodic --prefix test
    julia generate_ising_clusters_one_site.jl --help
"""

include("dependencies.jl")
include("functions/ClusterEnumeration.jl")

using ArgParse
using ProgressMeter

# Import lattice creation functions
function create_open_square_lattice(L::Int)
    """Create square lattice adjacency matrix with open boundary conditions."""
    N = L * L
    adj_matrix = zeros(Int, N, N)
    
    function coord_to_idx(i::Int, j::Int)
        return (i - 1) * L + j
    end
    
    for i in 1:L
        for j in 1:L
            current = coord_to_idx(i, j)
            
            # Right neighbor (only if not at right edge)
            if j < L
                right = coord_to_idx(i, j + 1)
                adj_matrix[current, right] = 1
                adj_matrix[right, current] = 1
            end
            
            # Down neighbor (only if not at bottom edge)
            if i < L
                down = coord_to_idx(i + 1, j)
                adj_matrix[current, down] = 1
                adj_matrix[down, current] = 1
            end
        end
    end
    
    return adj_matrix
end

function create_square_lattice(L::Int, boundary::String)
    """Create square lattice with specified boundary conditions."""
    if boundary == "periodic"
        return create_periodic_square_lattice(L)
    elseif boundary == "open"
        return create_open_square_lattice(L)
    else
        error("Unknown boundary condition: $boundary. Use 'periodic' or 'open'.")
    end
end

# Helper function to estimate completion time
function estimate_completion_time_one_site(L::Int, max_weight::Int)
    """Provide rough time estimates based on lattice size and weight for single site."""
    # Much faster than all-sites version
    if L <= 6
        return "< 5 seconds"
    elseif L <= 8
        return "< 30 seconds"
    elseif L <= 10
        return "1-5 minutes"
    elseif L <= 12
        return "5-15 minutes"
    else
        return "15+ minutes (consider smaller parameters)"
    end
end

# Distance and translation functions for square lattices
function site_to_coords(site::Int, L::Int)
    """Convert 1-indexed site to (i, j) coordinates."""
    i = div(site - 1, L) + 1
    j = mod(site - 1, L) + 1
    return (i, j)
end

function coords_to_site(i::Int, j::Int, L::Int)
    """Convert (i, j) coordinates to 1-indexed site."""
    return (i - 1) * L + j
end

function manhattan_distance_periodic(site1::Int, site2::Int, L::Int)
    """Calculate Manhattan distance between two sites on periodic lattice."""
    i1, j1 = site_to_coords(site1, L)
    i2, j2 = site_to_coords(site2, L)
    
    # Periodic distance in each dimension
    di = min(abs(i2 - i1), L - abs(i2 - i1))
    dj = min(abs(j2 - j1), L - abs(j2 - j1))
    
    return di + dj
end

function manhattan_distance_open(site1::Int, site2::Int, L::Int)
    """Calculate Manhattan distance between two sites on open lattice."""
    i1, j1 = site_to_coords(site1, L)
    i2, j2 = site_to_coords(site2, L)
    
    return abs(i2 - i1) + abs(j2 - j1)
end

function manhattan_distance(site1::Int, site2::Int, L::Int, boundary::String)
    """Calculate Manhattan distance with appropriate boundary conditions."""
    if boundary == "periodic"
        return manhattan_distance_periodic(site1, site2, L)
    else
        return manhattan_distance_open(site1, site2, L)
    end
end

function sites_within_distance(target_site::Int, max_distance::Int, L::Int, boundary::String)
    """Find all sites within Manhattan distance from target site."""
    nearby_sites = Int[]
    
    for site in 1:(L*L)
        if manhattan_distance(target_site, site, L, boundary) <= max_distance
            push!(nearby_sites, site)
        end
    end
    
    return nearby_sites
end

# Translation functions for periodic lattices
function translate_site(site::Int, di::Int, dj::Int, L::Int)
    """Translate a site by (di, dj) on an L×L periodic lattice."""
    i, j = site_to_coords(site, L)
    
    # Apply translation with periodic boundary conditions
    new_i = mod(i - 1 + di, L) + 1
    new_j = mod(j - 1 + dj, L) + 1
    
    return coords_to_site(new_i, new_j, L)
end

function translate_cluster(cluster::Cluster, di::Int, dj::Int, L::Int, all_loops::Vector{Loop})
    """Translate a cluster by (di, dj) and return the translated cluster."""
    # Create mapping from old loop to translated loop
    translated_loop_ids = Int[]
    translated_multiplicities = Dict{Int, Int}()
    
    # For each loop in the cluster, find its translation
    for (loop_id, multiplicity) in cluster.multiplicities
        original_loop = all_loops[loop_id]
        
        # Translate all vertices in the loop
        translated_vertices = [translate_site(v, di, dj, L) for v in original_loop.vertices]
        translated_edges = [(translate_site(u, di, dj, L), translate_site(v, di, dj, L)) for (u, v) in original_loop.edges]
        
        # Create canonical representation and find matching loop
        translated_loop = Loop(sort(translated_vertices), sort([(min(u,v), max(u,v)) for (u,v) in translated_edges]), original_loop.weight)
        translated_canonical = canonical_loop_representation(translated_loop)
        
        # Find the corresponding loop ID in all_loops
        translated_loop_id = -1
        for (i, loop) in enumerate(all_loops)
            if canonical_loop_representation(loop) == translated_canonical
                translated_loop_id = i
                break
            end
        end
        
        if translated_loop_id == -1
            error("Could not find translated loop - this should not happen!")
        end
        
        push!(translated_loop_ids, translated_loop_id)
        translated_multiplicities[translated_loop_id] = multiplicity
    end
    
    return Cluster(
        sort(unique(translated_loop_ids)),
        translated_multiplicities,
        cluster.weight,
        cluster.total_loops
    )
end

function canonical_cluster_under_translation(cluster::Cluster, L::Int, all_loops::Vector{Loop})
    """Find the canonical representative of a cluster under translation symmetry."""
    canonical_cluster = cluster
    canonical_signature = canonical_cluster_signature(cluster)
    
    # Try all possible translations
    for di in 0:(L-1)
        for dj in 0:(L-1)
            if di == 0 && dj == 0
                continue  # Skip identity translation
            end
            
            translated_cluster = translate_cluster(cluster, di, dj, L, all_loops)
            translated_signature = canonical_cluster_signature(translated_cluster)
            
            # Keep the lexicographically smallest signature
            if translated_signature < canonical_signature
                canonical_cluster = translated_cluster
                canonical_signature = translated_signature
            end
        end
    end
    
    return canonical_cluster
end

function remove_translational_redundancies(clusters::Vector{Cluster}, L::Int, all_loops::Vector{Loop})
    """Remove clusters that are related by translation, keeping only canonical representatives."""
    println("🔄 Removing translational redundancies...")
    flush(stdout)
    
    unique_clusters = Cluster[]
    seen_signatures = Set{Tuple}()
    
    progress = Progress(length(clusters), dt=0.1, desc="Removing redundancies: ", color=:orange, barlen=50)
    
    for cluster in clusters
        # Find the canonical representative under translation
        canonical_cluster = canonical_cluster_under_translation(cluster, L, all_loops)
        canonical_signature = canonical_cluster_signature(canonical_cluster)
        
        if !(canonical_signature in seen_signatures)
            push!(seen_signatures, canonical_signature)
            push!(unique_clusters, canonical_cluster)
        end
        
        next!(progress, showvalues = [("Unique", "$(length(unique_clusters))"), ("Total", "$(length(clusters))")])
    end
    
    finish!(progress)
    
    println("  Original clusters: $(length(clusters))")
    println("  Unique clusters (after translation): $(length(unique_clusters))")
    println("  Redundancy factor: $(round(length(clusters) / length(unique_clusters), digits=2))x")
    
    return unique_clusters
end

# Optimized loop enumeration for single site
function find_loops_within_distance(enumerator::ClusterEnumerator, target_site::Int, 
                                   max_weight::Int, max_distance::Int, L::Int, boundary::String)
    """Find all loops that contain vertices within max_distance of target_site."""
    
    println("  🎯 Optimizing: only considering vertices within distance $max_distance from site $target_site")
    
    # Find all sites within the distance constraint
    nearby_sites = sites_within_distance(target_site, max_distance, L, boundary)
    println("  Found $(length(nearby_sites)) sites within distance $max_distance (out of $(L*L) total)")
    println("  Reduction factor: $(round(L*L / length(nearby_sites), digits=2))x fewer sites to consider")
    
    # Create a set for fast lookup
    nearby_sites_set = Set(nearby_sites)
    
    all_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    # Progress bar for loop enumeration over nearby sites only
    progress = Progress(length(nearby_sites), dt=0.1, desc="Finding loops (optimized): ", color=:green, barlen=50)
    
    for vertex in nearby_sites
        # Find loops supported on this vertex
        vertex_loops = find_loops_supported_on_vertex(enumerator.loop_enumerator, vertex, max_weight)
        
        for loop in vertex_loops
            # Additional filter: only keep loops where ALL vertices are within our constraint
            # This ensures we don't miss any loops that could be part of relevant clusters
            if all(v in nearby_sites_set for v in loop.vertices)
                # Only keep loops where 'vertex' is the smallest vertex (canonical representative)
                if vertex == minimum(loop.vertices)
                    canonical = canonical_loop_representation(loop)
                    
                    if !(canonical in seen_loops)
                        push!(seen_loops, canonical)
                        new_loop = Loop(sort(loop.vertices), sort(loop.edges), loop.weight)
                        push!(all_loops, new_loop)
                    end
                end
            end
        end
        
        next!(progress, showvalues = [("Vertex", "$(findfirst(==(vertex), nearby_sites))/$(length(nearby_sites))"), ("Loops found", "$(length(all_loops))")])  
    end
    
    finish!(progress)
    
    return all_loops
end

# Data structure for single-site enumeration results
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
end

function enumerate_clusters_one_site_with_translation_removal(enumerator::ClusterEnumerator, site::Int, 
                                                           max_weight::Int, max_loop_weight::Int, L::Int, boundary::String)
    """Enumerate clusters for one site and remove translational redundancies."""
    
    println("🚀 Starting single-site cluster enumeration...")
    println("📋 Target site: $site")
    println("📋 Max cluster weight: $max_weight")
    println("📋 Max loop weight: $max_loop_weight")
    println("📋 Distance constraint: $max_weight (Manhattan distance)")
    println("📋 Boundary condition: $boundary")
    println()
    
    start_time = time()
    
    # Step 1: Enumerate loops within distance constraint (OPTIMIZATION)
    println("📊 Step 1: Enumerating loops within distance constraint...")
    println("   🚀 OPTIMIZATION: Only considering loops near target site")
    flush(stdout)
    
    # Use max_weight as the distance constraint since clusters can't extend beyond this
    max_distance = max_weight
    all_loops = find_loops_within_distance(enumerator, site, max_loop_weight, max_distance, L, boundary)
    println("Found $(length(all_loops)) loops within distance $max_distance")
    
    # Step 2: Build interaction graph
    println("\n🔗 Step 2: Building interaction graph...")
    flush(stdout)
    interaction_graph = build_interaction_graph_optimized(all_loops)
    
    # Step 3: Enumerate clusters for the target site
    println("\n🌐 Step 3: Enumerating clusters for site $site...")
    flush(stdout)
    
    # Find loops supported on target site
    supported_loop_ids = Int[]
    for (i, loop) in enumerate(all_loops)
        if site in loop.vertices
            push!(supported_loop_ids, i)
        end
    end
    println("  Found $(length(supported_loop_ids)) loops supported on site $site")
    
    if isempty(supported_loop_ids)
        println("  No loops found on site $site - returning empty results")
        return SingleSiteClusterData(
            Cluster[], all_loops, enumerator.adj_matrix, max_weight, L, site,
            0.0, 0.0, string(now()), "", 0, 0
        )
    end
    
    # DFS cluster enumeration
    site_clusters = dfs_enumerate_clusters_from_supported(all_loops, supported_loop_ids, 
                                                        max_weight, interaction_graph)
    
    # Remove standard redundancies first
    println("\n🧹 Step 4: Removing standard redundancies...")
    unique_clusters = remove_cluster_redundancies(site_clusters)
    clusters_before_translation = length(unique_clusters)
    
    enum_time = time() - start_time
    
    # Step 5: Remove translational redundancies (only for periodic boundary conditions)
    println("\n🔄 Step 5: Removing translational redundancies...")
    translation_start_time = time()
    
    clusters_after_translation = remove_translational_redundancies(unique_clusters, L, all_loops)
    clusters_after_count = length(clusters_after_translation)
    
    translation_time = time() - translation_start_time
    
    return clusters_before_translation, clusters_after_translation, clusters_after_count, 
           enum_time, translation_time, all_loops
end

function save_single_site_cluster_data(data::SingleSiteClusterData, prefix::String = "")
    """Save single-site cluster enumeration data to file."""
    
    # Create directory if it doesn't exist
    save_dir = "saved_clusters"
    if !isdir(save_dir)
        mkpath(save_dir)
        println("📁 Created directory: $save_dir")
    end
    
    # Generate filename
    timestamp = replace(data.timestamp, ":" => "-", "." => "-")
    base_name = "single_site_clusters_L$(data.lattice_size)_site$(data.site)_w$(data.max_weight)_$(data.boundary_condition)_$(timestamp)"
    if !isempty(prefix)
        base_name = "$(prefix)_$(base_name)"
    end
    
    filepath = joinpath(save_dir, "$(base_name).jld2")
    
    # Save using Julia's serialization
    println("💾 Saving single-site enumeration data to: $(filepath)")
    
    # Create summary for quick inspection
    summary = Dict(
        "lattice_size" => data.lattice_size,
        "site" => data.site,
        "max_weight" => data.max_weight,
        "boundary_condition" => data.boundary_condition,
        "total_loops" => length(data.all_loops),
        "clusters_before_translation_removal" => data.clusters_before_translation_removal,
        "clusters_after_translation_removal" => data.clusters_after_translation_removal,
        "redundancy_factor" => data.clusters_before_translation_removal / max(1, data.clusters_after_translation_removal),
        "enumeration_time_seconds" => data.enumeration_time,
        "translation_removal_time_seconds" => data.translation_removal_time,
        "total_time_seconds" => data.enumeration_time + data.translation_removal_time,
        "timestamp" => data.timestamp
    )
    
    # Save both data and summary
    save_data = Dict(
        "data" => data,
        "summary" => summary
    )
    
    open(filepath, "w") do io
        serialize(io, save_data)
    end
    
    println("✅ Saved successfully!")
    println("   Summary: $summary")
    
    return filepath
end

function parse_commandline()
    s = ArgParseSettings(
        description = "Generate connected clusters for single site on square lattice Ising model",
        epilog = "Example: julia generate_ising_clusters_one_site.jl --size 6 --weight 4 --site 1 --boundary periodic"
    )

    @add_arg_table! s begin
        "--size", "-L"
            help = "Linear lattice size (L×L lattice)"
            arg_type = Int
            default = 6
            
        "--weight", "-w"
            help = "Maximum cluster weight to enumerate"
            arg_type = Int
            default = 4
            
        "--site", "-s"
            help = "Target site (1-indexed)"
            arg_type = Int
            default = 1
            
        "--loop-weight"
            help = "Maximum individual loop weight (defaults to cluster weight)"
            arg_type = Int
            default = nothing
            
        "--boundary", "-b"
            help = "Boundary conditions: 'periodic' or 'open'"
            arg_type = String
            default = "periodic"
            
        "--prefix", "-p"
            help = "Optional prefix for saved files"
            arg_type = String
            default = ""
            
        "--list"
            help = "List existing saved cluster files and exit"
            action = :store_true
            
        "--analyze"
            help = "Analyze a saved cluster file"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    
    println("🎯 Single-Site Ising Cluster Generator")
    println("="^50)
    
    # Handle special actions
    if args["list"]
        list_saved_cluster_files()
        return
    end
    
    if !isempty(args["analyze"])
        filepath = args["analyze"]
        if !isfile(filepath)
            # Try looking in saved_clusters directory
            filepath = joinpath("saved_clusters", filepath)
            if !isfile(filepath)
                println("❌ File not found: $(args["analyze"])")
                return
            end
        end
        
        println("📊 Analyzing: $filepath")
        loaded_data = open(filepath, "r") do io
            deserialize(io)
        end
        
        if haskey(loaded_data, "summary")
            summary = loaded_data["summary"]
            println("✅ Loaded successfully!")
            println("   Summary: $summary")
        else
            println("❌ Invalid file format")
        end
        return
    end
    
    # Extract parameters
    L = args["size"]
    max_weight = args["weight"]
    site = args["site"]
    boundary = args["boundary"]
    prefix = args["prefix"]
    max_loop_weight = something(args["loop-weight"], max_weight)
    
    # Validate parameters
    if L < 2
        println("❌ Error: Lattice size must be at least 2")
        return
    end
    
    if max_weight < 1
        println("❌ Error: Max weight must be at least 1")
        return
    end
    
    if site < 1 || site > L*L
        println("❌ Error: Site must be between 1 and $(L*L)")
        return
    end
    
    if !(boundary in ["periodic", "open"])
        println("❌ Error: Boundary must be 'periodic' or 'open'")
        return
    end
    
    # Display parameters with time estimate
    println("📋 Parameters:")
    println("  Lattice size: $(L)×$(L) (L^2 = $(L*L) sites)")
    println("  Target site: $site")
    println("  Boundary conditions: $boundary")
    println("  Max cluster weight: $max_weight")
    println("  Max loop weight: $max_loop_weight")
    if !isempty(prefix)
        println("  File prefix: '$prefix'")
    end
    
    # Add time estimate
    time_estimate = estimate_completion_time_one_site(L, max_weight)
    println("  ⏱️  Estimated time: $time_estimate")
    
    if boundary == "periodic"
        println("  🔄 Translation symmetries will be removed")
    else
        println("  ⚠️  Open boundaries: translation removal not applicable")
    end
    
    # Show optimization info
    total_sites = L * L
    max_distance = max_weight
    nearby_sites_count = length(sites_within_distance(site, max_distance, L, boundary))
    optimization_factor = round(total_sites / nearby_sites_count, digits=2)
    
    println("  🎯 OPTIMIZATION INFO:")
    println("    Distance constraint: $max_distance (based on max weight)")
    println("    Sites within constraint: $nearby_sites_count / $total_sites")
    println("    Expected speedup: ~$(optimization_factor)x faster loop enumeration")
    println()
    
    # Create lattice with progress
    println("🏗️  Creating $(L)×$(L) square lattice with $boundary boundary conditions...")
    lattice_progress = Progress(3, dt=0.5, desc="Setting up lattice: ", color=:blue)
    
    next!(lattice_progress, showvalues = [("Step", "Creating adjacency matrix")])
    adj_matrix = create_square_lattice(L, boundary)
    
    next!(lattice_progress, showvalues = [("Step", "Initializing enumerator")])
    enumerator = ClusterEnumerator(adj_matrix)
    
    next!(lattice_progress, showvalues = [("Step", "Lattice setup complete")])
    finish!(lattice_progress)
    
    # Add boundary condition to prefix
    full_prefix = isempty(prefix) ? boundary : "$(prefix)_$(boundary)"
    
    start_time = time()
    
    try
        # Enumerate clusters for the single site
        clusters_before_count, clusters_after_translation, clusters_after_count, 
        enum_time, translation_time, all_loops = enumerate_clusters_one_site_with_translation_removal(
            enumerator, site, max_weight, max_loop_weight, L, boundary)
        
        total_time = time() - start_time
        
        # Create data structure
        data = SingleSiteClusterData(
            clusters_after_translation,
            all_loops,
            adj_matrix,
            max_weight,
            L,
            site,
            enum_time,
            translation_time,
            string(now()),
            boundary,
            clusters_before_count,
            clusters_after_count
        )
        
        # Save results
        filepath = save_single_site_cluster_data(data, full_prefix)
        
        # Success summary
        println("\n🎉 SUCCESS!")
        println("="^40)
        println("✅ Single-site enumeration completed in $(round(total_time, digits=2)) seconds")
        println("✅ Clusters before translation removal: $clusters_before_count")
        println("✅ Clusters after translation removal: $clusters_after_count")
        if clusters_before_count > 0
            reduction_factor = round(clusters_before_count / max(1, clusters_after_count), digits=2)
            println("✅ Redundancy reduction factor: $(reduction_factor)x")
        end
        println("✅ Data saved successfully!")
        println("📁 Saved to: $(basename(filepath))")
        
        # Quick analysis
        println("\n📊 Quick Analysis:")
        by_weight = Dict{Int, Int}()
        for cluster in clusters_after_translation
            by_weight[cluster.weight] = get(by_weight, cluster.weight, 0) + 1
        end
        
        for weight in sort(collect(keys(by_weight)))
            count = by_weight[weight]
            println("  Weight $weight: $count unique clusters")
        end
        
    catch e
        println("\n❌ ERROR during enumeration:")
        println("   $e")
        
        if isa(e, InterruptException)
            println("   (Interrupted by user)")
        else
            println("   Check parameters and try again")
        end
        
        return
    end
    
    println("\n💡 Next steps:")
    println("  • Use --list to see all saved files")
    println("  • Use --analyze <filename> to analyze results")
    println("  • Compare with all-sites enumeration using generate_ising_clusters.jl")
end

function show_examples()
    """Show usage examples."""
    println("📚 EXAMPLES:")
    println()
    println("1. Generate clusters for site 1 on 4×4 periodic lattice, max weight 3:")
    println("   julia generate_ising_clusters_one_site.jl --size 4 --site 1 --weight 3 --boundary periodic")
    println()
    println("2. Generate clusters for site 10 on 6×6 open lattice with custom prefix:")
    println("   julia generate_ising_clusters_one_site.jl -L 6 -s 10 -w 5 -b open -p my_test")
    println()
    println("3. List all saved cluster files:")
    println("   julia generate_ising_clusters_one_site.jl --list")
    println()
    println("4. Analyze a saved file:")
    println("   julia generate_ising_clusters_one_site.jl --analyze saved_clusters/single_site_*.jld2")
    println()
end

# Handle help and examples
if length(ARGS) > 0 && (ARGS[1] == "--help" || ARGS[1] == "-h")
    # ArgParse will handle the basic help, but we can add examples
    try
        parse_commandline()
    catch e
        if isa(e, ArgParseError)
            # This will show the help
            rethrow(e)
        end
    end
    println()
    show_examples()
    exit(0)
elseif length(ARGS) > 0 && ARGS[1] == "--examples"
    show_examples()
    exit(0)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end