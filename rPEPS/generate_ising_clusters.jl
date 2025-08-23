#!/usr/bin/env julia

"""
generate_ising_clusters.jl

Generate and save connected cluster enumerations for square lattice Ising models.
Supports both periodic and open boundary conditions.

Usage:
    julia generate_ising_clusters.jl --size 6 --weight 4 --boundary periodic --prefix test
    julia generate_ising_clusters.jl --help
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")
include("../functions/loop_enumeration_square.jl")

using ArgParse
using ProgressMeter

# Helper function to estimate completion time
function estimate_completion_time(L::Int, max_weight::Int)
    """Provide rough time estimates based on lattice size and weight."""
    # Updated estimates for the new optimized approach
    if L <= 4
        return "< 1 second"
    elseif L <= 6
        return "< 10 seconds"
    elseif L <= 8
        return "10-60 seconds"
    elseif L <= 10
        return "1-5 minutes"
    else
        return "5+ minutes (consider smaller parameters)"
    end
end

# Center site calculation
function calculate_center_site(L::Int)
    """Calculate the center site index for an L×L lattice."""
    center_i = (L + 1) ÷ 2
    center_j = (L + 1) ÷ 2
    return (center_i - 1) * L + center_j
end

# Translation functions for loops
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

function translate_site(site::Int, di::Int, dj::Int, L::Int)
    """Translate a site by (di, dj) on an L×L periodic lattice."""
    i, j = site_to_coords(site, L)
    
    # Apply translation with periodic boundary conditions
    new_i = mod(i - 1 + di, L) + 1
    new_j = mod(j - 1 + dj, L) + 1
    
    return coords_to_site(new_i, new_j, L)
end

function translate_loop(loop::Loop, di::Int, dj::Int, L::Int)
    """Translate a loop by (di, dj) and return the translated loop."""
    # Translate all vertices in the loop
    translated_vertices = [translate_site(v, di, dj, L) for v in loop.vertices]
    translated_edges = [(translate_site(u, di, dj, L), translate_site(v, di, dj, L)) for (u, v) in loop.edges]
    
    # Create canonical representation
    translated_loop = Loop(sort(translated_vertices), sort([(min(u,v), max(u,v)) for (u,v) in translated_edges]), loop.weight)
    return translated_loop
end

function is_loop_within_bounds(loop::Loop, L::Int)
    """Check if all vertices of a loop are within [1, L²] bounds."""
    return all(1 <= v <= L*L for v in loop.vertices)
end

function generate_all_loops_pbc(center_loops::Vector{Loop}, L::Int)
    """Generate all loops by translating center loops for periodic boundary conditions."""
    all_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    println("  🔄 Generating all loops via translation (PBC)...")
    total_translations = L * L
    progress = Progress(total_translations, dt=0.1, desc="Translating loops: ", color=:green, barlen=50)
    
    for di in 0:(L-1)
        for dj in 0:(L-1)
            for loop in center_loops
                translated_loop = translate_loop(loop, di, dj, L)
                canonical = canonical_loop_representation(translated_loop)
                
                if !(canonical in seen_loops)
                    push!(seen_loops, canonical)
                    push!(all_loops, translated_loop)
                end
            end
            next!(progress, showvalues = [("Translation", "($di,$dj)"), ("Unique loops", "$(length(all_loops))")])
        end
    end
    
    finish!(progress)
    println("  Generated $(length(all_loops)) unique loops from $(length(center_loops)) center loops")
    
    return all_loops
end

function generate_all_loops_obc(center_loops::Vector{Loop}, L::Int)
    """Generate all loops by translating center loops for open boundary conditions."""
    all_loops = Loop[]
    seen_loops = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    println("  🔄 Generating all loops via translation (OBC)...")
    total_translations = L * L
    progress = Progress(total_translations, dt=0.1, desc="Translating loops: ", color=:green, barlen=50)
    
    discarded_count = 0
    
    for di in 0:(L-1)
        for dj in 0:(L-1)
            for loop in center_loops
                translated_loop = translate_loop(loop, di, dj, L)
                
                # Check if translated loop is within bounds
                if is_loop_within_bounds(translated_loop, L)
                    canonical = canonical_loop_representation(translated_loop)
                    
                    if !(canonical in seen_loops)
                        push!(seen_loops, canonical)
                        push!(all_loops, translated_loop)
                    end
                else
                    discarded_count += 1
                end
            end
            next!(progress, showvalues = [("Translation", "($di,$dj)"), ("Unique loops", "$(length(all_loops))"), ("Discarded", "$discarded_count")])
        end
    end
    
    finish!(progress)
    println("  Generated $(length(all_loops)) unique loops from $(length(center_loops)) center loops")
    println("  Discarded $discarded_count out-of-bounds translations")
    
    return all_loops
end


function create_open_square_lattice(L::Int)
    """Create square lattice adjacency matrix with open boundary conditions."""
    N = L * L
    adj_matrix = zeros(Int, N, N)
    
    function coord_to_idx(i::Int, j::Int)
        return (i - 1) * L + j
    end
    
    # Add progress bar for lattice construction if it's a larger lattice
    if L >= 6
        progress = Progress(L, dt=0.1, desc="Building lattice: ", color=:blue)
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
            next!(progress, showvalues = [("Row", "$i/$L")])
        end
        finish!(progress)
    else
        # For small lattices, just build without progress bar
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

function parse_commandline()
    s = ArgParseSettings(
        description = "Generate connected clusters for square lattice Ising model",
        epilog = "Example: julia generate_ising_clusters.jl --size 6 --weight 4 --boundary periodic"
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

function create_canonical_cluster_signature(cluster::Cluster, all_loops::Vector{Loop}, L::Int = 11)
    """
    Create cluster signature using canonical loop representations instead of loop IDs.
    This fixes the deduplication issue where equivalent loops have different IDs.
    """
    canonical_loop_sigs = []
    for loop_id in cluster.loop_ids
        loop = all_loops[loop_id]
        # Create proper coordinate-based canonical form
        canonical_loop = create_coordinate_canonical_form(loop, L)
        multiplicity = cluster.multiplicities[loop_id]
        push!(canonical_loop_sigs, (canonical_loop, multiplicity))
    end
    
    # Sort by canonical loop representation for consistent ordering
    sort!(canonical_loop_sigs)
    
    return (tuple(canonical_loop_sigs...), cluster.weight)
end

function create_coordinate_canonical_form(loop::Loop, L::Int = 11)
    """
    Create a true canonical form based on coordinate geometry.
    This handles periodic boundary conditions properly for all pattern types.
    """
    # Convert vertices to coordinates
    function site_to_coords(site::Int)
        i = div(site - 1, L) + 1
        j = mod(site - 1, L) + 1
        return (i, j)
    end
    
    coords = [site_to_coords(v) for v in sort(loop.vertices)]
    
    # Check for periodic wrapping by examining coordinate spans
    i_values = [coord[1] for coord in coords]
    j_values = [coord[2] for coord in coords]
    
    i_span = maximum(i_values) - minimum(i_values)
    j_span = maximum(j_values) - minimum(j_values)
    
    # If span > L/2, the pattern likely wraps around periodic boundaries
    # We need to "unwrap" it to get the true geometric pattern
    if i_span > L/2 || j_span > L/2
        # Unwrap coordinates by shifting wrapped vertices
        unwrapped_coords = []
        
        for coord in coords
            i, j = coord
            # Unwrap i-coordinate if needed
            if i_span > L/2
                # Find vertices that are "far" from the minimum i
                min_i = minimum(i_values)
                if i - min_i > L/2
                    i = i - L  # Shift back to represent wrapping
                end
            end
            
            # Unwrap j-coordinate if needed  
            if j_span > L/2
                # Find vertices that are "far" from the minimum j
                min_j = minimum(j_values)
                if j - min_j > L/2
                    j = j - L  # Shift back to represent wrapping
                end
            end
            
            push!(unwrapped_coords, (i, j))
        end
        
        coords = unwrapped_coords
    end
    
    # Normalize coordinates to start from (0, 0)
    min_i = minimum(coord[1] for coord in coords)
    min_j = minimum(coord[2] for coord in coords)
    
    normalized_coords = [(coord[1] - min_i, coord[2] - min_j) for coord in coords]
    
    # Sort for canonical ordering
    normalized_coords = sort(normalized_coords)
    
    return (normalized_coords, loop.weight)
end

function apply_global_cluster_deduplication(clusters_by_site::Dict{Int, Vector{Cluster}}, all_loops::Vector{Loop}, L::Int)
    """
    Remove cluster duplicates across all sites using canonical signatures.
    This handles cases where the same cluster pattern appears multiple times
    due to loop enumeration creating equivalent loops with different IDs.
    """
    
    println("  🔍 Analyzing cluster duplicates...")
    
    # Collect all clusters with their sites
    all_clusters_with_sites = []
    original_total = 0
    
    for (site, clusters) in clusters_by_site
        original_total += length(clusters)
        for cluster in clusters
            push!(all_clusters_with_sites, (site, cluster))
        end
    end
    
    println("  📊 Processing $(length(all_clusters_with_sites)) total clusters from $original_total site assignments")
    
    # Group by canonical cluster signature (using loop canonical forms)
    canonical_groups = Dict{Tuple, Vector{Tuple{Int, Cluster}}}()
    
    progress = Progress(length(all_clusters_with_sites), dt=0.1, desc="  Creating canonical signatures: ", color=:yellow, barlen=40)
    
    for (site, cluster) in all_clusters_with_sites
        # Create signature using canonical loop representations
        canonical_sig = create_canonical_cluster_signature(cluster, all_loops, L)
        
        if !haskey(canonical_groups, canonical_sig)
            canonical_groups[canonical_sig] = []
        end
        push!(canonical_groups[canonical_sig], (site, cluster))
        
        next!(progress)
    end
    
    finish!(progress)
    
    println("  📊 Found $(length(canonical_groups)) unique canonical cluster forms")
    
    # For each canonical group, place one copy on each site where it appears
    deduplicated_clusters_by_site = Dict{Int, Vector{Cluster}}()
    for site in 1:(L*L)
        deduplicated_clusters_by_site[site] = []
    end
    
    dedup_progress = Progress(length(canonical_groups), dt=0.1, desc="  Deduplicating: ", color=:yellow, barlen=40)
    
    total_duplicates_removed = 0
    
    for (canonical_sig, site_cluster_pairs) in canonical_groups
        # Group by site to see how many copies appear on each site
        clusters_by_site_for_this_canonical = Dict{Int, Vector{Cluster}}()
        for (site, cluster) in site_cluster_pairs
            if !haskey(clusters_by_site_for_this_canonical, site)
                clusters_by_site_for_this_canonical[site] = []
            end
            push!(clusters_by_site_for_this_canonical[site], cluster)
        end
        
        # For each site, keep only one copy (remove duplicates)
        for (site, clusters_on_site) in clusters_by_site_for_this_canonical
            # Keep only the first cluster (they're all equivalent)
            representative_cluster = clusters_on_site[1]
            push!(deduplicated_clusters_by_site[site], representative_cluster)
            
            # Count how many duplicates we removed from this site
            duplicates_removed_from_site = length(clusters_on_site) - 1
            total_duplicates_removed += duplicates_removed_from_site
        end
        
        next!(dedup_progress)
    end
    
    finish!(dedup_progress)
    
    # Calculate final statistics
    final_total = sum(length(clusters) for clusters in values(deduplicated_clusters_by_site))
    
    println("  ✅ Deduplication complete:")
    println("    Original clusters: $original_total")
    println("    Final clusters: $final_total") 
    println("    Duplicates removed: $(original_total - final_total)")
    println("    Unique canonical forms: $(length(canonical_groups))")
    
    return deduplicated_clusters_by_site
end

function main()
    args = parse_commandline()
    
    println("🎯 Ising Cluster Generator")
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
        data = load_cluster_enumeration(filepath)
        analyze_saved_enumeration(data)
        return
    end
    
    # Extract parameters
    L = args["size"]
    max_weight = args["weight"]
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
    
    if !(boundary in ["periodic", "open"])
        println("❌ Error: Boundary must be 'periodic' or 'open'")
        return
    end
    
    # Display parameters with time estimate
    println("📋 Parameters:")
    println("  Lattice size: $(L)×$(L) (L^2 = $(L*L) sites)")
    println("  Boundary conditions: $boundary")
    println("  Max cluster weight: $max_weight")
    println("  Max loop weight: $max_loop_weight")
    if !isempty(prefix)
        println("  File prefix: '$prefix'")
    end
    
    # Add time estimate
    time_estimate = estimate_completion_time(L, max_weight)
    println("  ⏱️  Estimated time: $time_estimate")
    println()
    
    # Note: Lattice creation is now handled by SquareLoopEnumerator during enumeration
    println("🏗️  $(L)×$(L) square lattice will be created during enumeration with $boundary boundary conditions...")
    
    # Add boundary condition to prefix
    full_prefix = isempty(prefix) ? boundary : "$(prefix)_$(boundary)"
    
    # Generate clusters using the new optimized approach
    println("\n🚀 Starting optimized cluster enumeration...")
    println("📊 This will involve 4 main steps:")
    println("   1️⃣  Center site loop generation")
    println("   2️⃣  Loop translation and expansion")
    println("   3️⃣  Cluster enumeration per site")
    println("   4️⃣  Data saving and analysis")
    println()
    
    start_time = time()
    
    try
        # Step 1: Generate loops supported on center site
        center_site = calculate_center_site(L)
        center_coords = site_to_coords(center_site, L)
        println("📍 Step 1: Generating loops on center site $center_site (coordinates $center_coords)...")
        
        # Use SquareLoopEnumerator for optimized performance
        enumerator_square = SquareLoopEnumerator(L; periodic=(boundary == "periodic"))
        center_loops = find_loops_supported_on_vertex_square(enumerator_square, center_site, max_loop_weight)
        
        println("  Found $(length(center_loops)) loops supported on center site")
        
        # Step 2: Generate all loops via translation
        println("\n🔄 Step 2: Generating all loops via translation...")
        if boundary == "periodic"
            all_loops = generate_all_loops_pbc(center_loops, L)
        else
            all_loops = generate_all_loops_obc(center_loops, L)
        end
        
        # Step 3: Build interaction graph and enumerate clusters
        println("\n🔗 Step 3: Building interaction graph and enumerating clusters...")
        interaction_graph = build_interaction_graph_optimized(all_loops)
        
        # Enumerate clusters for each site (following one_site_old.jl pattern)
        clusters_by_site = Dict{Int, Vector{Cluster}}()
        n_sites = L * L
        
        println("  Processing $(n_sites) sites...")
        site_progress = Progress(n_sites, dt=0.1, desc="Sites processed: ", color=:cyan, barlen=50)
        
        for site in 1:n_sites
            # Find loops supported on this site
            supported_loop_ids = Int[]
            for (i, loop) in enumerate(all_loops)
                if site in loop.vertices
                    push!(supported_loop_ids, i)
                end
            end
            
            if !isempty(supported_loop_ids)
                # Enumerate clusters for this site
                site_clusters = dfs_enumerate_clusters_from_supported(all_loops, supported_loop_ids, 
                                                                    max_weight, interaction_graph)
                # Remove redundancies
                unique_clusters = remove_cluster_redundancies(site_clusters)
                clusters_by_site[site] = unique_clusters
            else
                clusters_by_site[site] = Cluster[]
            end
            
            next!(site_progress)
        end
        
        finish!(site_progress)
        
        # Step 3.5: Global cluster deduplication across all sites
        println("\n🔧 Step 3.5: Global cluster deduplication...")
        clusters_by_site = apply_global_cluster_deduplication(clusters_by_site, all_loops, L)
        
        # Step 4: Create data structure and save
        println("\n💾 Step 4: Creating data structure and saving...")
        enumeration_time = time() - start_time
        
        data = ClusterEnumerationData(
            clusters_by_site,
            all_loops,
            enumerator_square.base_enumerator.adj_matrix,
            max_weight,
            L,
            enumeration_time,
            string(now()),
            n_sites
        )
        
        # Save results
        save_cluster_enumeration(data, full_prefix)
        
        println("\n✅ Core enumeration completed!")
        
        end_time = time()
        total_time = end_time - start_time
        
        # Success summary with progress
        println("⏳ Calculating final statistics...")
        stats_progress = Progress(3, dt=0.2, desc="Final analysis: ", color=:cyan)
        
        next!(stats_progress, showvalues = [("Task", "Counting clusters")])
        total_clusters = sum(length(clusters) for clusters in values(data.clusters_by_site))
        
        next!(stats_progress, showvalues = [("Task", "Computing averages")])
        # Give time for the progress bar to show
        sleep(0.1)
        
        next!(stats_progress, showvalues = [("Task", "Preparing output")])
        finish!(stats_progress)
        
        println("\n🎉 SUCCESS!")
        println("="^40)
        println("✅ Enumeration completed in $(round(total_time, digits=2)) seconds")
        println("✅ Total clusters found: $total_clusters")
        println("✅ Average per site: $(round(total_clusters / data.total_sites, digits=2))")
        println("✅ Data saved successfully!")
        
        # Show file location with progress
        file_progress = Progress(2, dt=0.2, desc="Locating files: ", color=:magenta)
        
        next!(file_progress, showvalues = [("Task", "Scanning directory")])
        save_dir = "saved_clusters"
        files = filter(f -> startswith(f, full_prefix) && endswith(f, ".jld2"), 
                      readdir(save_dir; join=false))
        
        next!(file_progress, showvalues = [("Task", "Finding latest file")])
        finish!(file_progress)
        
        if !isempty(files)
            latest_file = sort(files)[end]  # Get most recent
            println("📁 Saved to: $save_dir/$latest_file")
        end
        
        # Quick analysis with progress
        println("\n📊 Quick Analysis:")
        analysis_progress = Progress(2, dt=0.3, desc="Analyzing results: ", color=:yellow)
        
        next!(analysis_progress, showvalues = [("Task", "Grouping by weight")])
        by_weight = Dict{Int, Int}()
        for clusters in values(data.clusters_by_site)
            for cluster in clusters
                by_weight[cluster.weight] = get(by_weight, cluster.weight, 0) + 1
            end
        end
        
        next!(analysis_progress, showvalues = [("Task", "Generating summary")])
        finish!(analysis_progress)
        
        for weight in sort(collect(keys(by_weight)))
            count = by_weight[weight]
            println("  Weight $weight: $count clusters")
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
    println("  • Load data in Julia with: load_cluster_enumeration(\"path/to/file.jld2\")")
end

function show_examples()
    """Show usage examples."""
    println("📚 EXAMPLES:")
    println()
    println("1. Generate clusters for 4×4 periodic lattice, max weight 3:")
    println("   julia generate_ising_clusters.jl --size 4 --weight 3 --boundary periodic")
    println()
    println("2. Generate clusters for 6×6 open lattice with custom prefix:")
    println("   julia generate_ising_clusters.jl -L 6 -w 5 -b open -p my_test")
    println()
    println("3. List all saved cluster files:")
    println("   julia generate_ising_clusters.jl --list")
    println()
    println("4. Analyze a saved file:")
    println("   julia generate_ising_clusters.jl --analyze saved_clusters/periodic_clusters_L6_w4_*.jld2")
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