#!/usr/bin/env julia

"""
generate_ising_clusters_one_site.jl

Generate and save connected cluster enumerations for a single site on square lattice Ising models.
Only supports periodic boundary conditions with translation symmetry deduplication.
Uses the same functions as generate_ising_clusters.jl for consistency.

Usage:
    julia generate_ising_clusters_one_site.jl --size 6 --weight 4 --site 1 --prefix test
    julia generate_ising_clusters_one_site.jl --help
"""

include("dependencies.jl")
include("functions/ClusterEnumeration.jl")
include("functions/loop_enumeration_square.jl")

using ArgParse
using ProgressMeter

# Import functions from generate_ising_clusters.jl
function estimate_completion_time(L::Int, max_weight::Int)
    """Provide rough time estimates based on lattice size and weight."""
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

function calculate_center_site(L::Int)
    """Calculate the center site index for an L×L lattice."""
    center_i = (L + 1) ÷ 2
    center_j = (L + 1) ÷ 2
    return (center_i - 1) * L + center_j
end

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
    if i_span > L/2 || j_span > L/2
        # Unwrap coordinates by shifting wrapped vertices
        unwrapped_coords = []
        
        for coord in coords
            i, j = coord
            # Unwrap i-coordinate if needed
            if i_span > L/2
                min_i = minimum(i_values)
                if i - min_i > L/2
                    i = i - L  # Shift back to represent wrapping
                end
            end
            
            # Unwrap j-coordinate if needed  
            if j_span > L/2
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
    normalized_coords = sort(normalized_coords)
    
    return (normalized_coords, loop.weight)
end

function create_translation_aware_cluster_signature(cluster::Cluster, all_loops::Vector{Loop}, L::Int)
    """
    Create cluster signature using translation-aware canonical loop forms.
    This ensures clusters related by translation are identified as equivalent.
    """
    canonical_loop_sigs = []
    for loop_id in cluster.loop_ids
        loop = all_loops[loop_id]
        # Create coordinate-based canonical form for translation equivalence
        canonical_loop = create_coordinate_canonical_form(loop, L)
        multiplicity = cluster.multiplicities[loop_id]
        push!(canonical_loop_sigs, (canonical_loop, multiplicity))
    end
    
    # Sort by canonical loop representation for consistent ordering
    sort!(canonical_loop_sigs)
    
    return (tuple(canonical_loop_sigs...), cluster.weight)
end

function apply_translation_aware_deduplication(clusters::Vector{Cluster}, all_loops::Vector{Loop}, L::Int)
    """
    Remove clusters that are equivalent under translation symmetry.
    Keep only one representative from each translation equivalence class.
    """
    
    println("  🔍 Applying translation-aware deduplication...")
    
    # Group clusters by translation-aware canonical signature
    canonical_groups = Dict{Tuple, Vector{Cluster}}()
    
    progress = Progress(length(clusters), dt=0.1, desc="  Creating canonical signatures: ", color=:yellow, barlen=40)
    
    for cluster in clusters
        # Create translation-aware signature
        canonical_sig = create_translation_aware_cluster_signature(cluster, all_loops, L)
        
        if !haskey(canonical_groups, canonical_sig)
            canonical_groups[canonical_sig] = []
        end
        push!(canonical_groups[canonical_sig], cluster)
        
        next!(progress)
    end
    
    finish!(progress)
    
    println("  📊 Found $(length(canonical_groups)) unique canonical cluster forms")
    
    # Keep only one representative from each equivalence class
    deduplicated_clusters = Cluster[]
    
    dedup_progress = Progress(length(canonical_groups), dt=0.1, desc="  Selecting representatives: ", color=:yellow, barlen=40)
    
    total_duplicates_removed = 0
    
    for (canonical_sig, equivalent_clusters) in canonical_groups
        # Keep the first cluster (they're all equivalent under translation)
        representative_cluster = equivalent_clusters[1]
        push!(deduplicated_clusters, representative_cluster)
        
        # Count duplicates removed
        duplicates_removed = length(equivalent_clusters) - 1
        total_duplicates_removed += duplicates_removed
        
        next!(dedup_progress)
    end
    
    finish!(dedup_progress)
    
    println("  ✅ Translation-aware deduplication complete:")
    println("    Original clusters: $(length(clusters))")
    println("    Final clusters: $(length(deduplicated_clusters))") 
    println("    Duplicates removed: $(total_duplicates_removed)")
    println("    Unique canonical forms: $(length(canonical_groups))")
    println("    Reduction factor: $(round(length(clusters) / length(deduplicated_clusters), digits=2))x")
    
    return deduplicated_clusters
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
    canonical_forms_count::Int
end

function parse_commandline()
    s = ArgParseSettings(
        description = "Generate connected clusters for single site on square lattice Ising model (periodic BC only)",
        epilog = "Example: julia generate_ising_clusters_one_site.jl --size 6 --weight 4 --site 1"
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
            help = "Target site (1-indexed, default is center site)"
            arg_type = Int
            default = nothing
            
        "--loop-weight"
            help = "Maximum individual loop weight (defaults to cluster weight)"
            arg_type = Int
            default = nothing
            
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
    base_name = "single_site_clusters_L$(data.lattice_size)_site$(data.site)_w$(data.max_weight)_periodic_$(timestamp)"
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
        "canonical_forms_count" => data.canonical_forms_count,
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

function main()
    args = parse_commandline()
    
    println("🎯 Single-Site Ising Cluster Generator (Translation-Aware)")
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
    site = something(args["site"], calculate_center_site(L))
    prefix = args["prefix"]
    max_loop_weight = something(args["loop-weight"], max_weight)
    boundary = "periodic"  # Only periodic boundary conditions supported
    
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
    
    # Display parameters with time estimate
    println("📋 Parameters:")
    println("  Lattice size: $(L)×$(L) (L^2 = $(L*L) sites)")
    println("  Target site: $site (coordinates $(site_to_coords(site, L)))")
    println("  Boundary conditions: $boundary (only supported option)")
    println("  Max cluster weight: $max_weight")
    println("  Max loop weight: $max_loop_weight")
    if !isempty(prefix)
        println("  File prefix: '$prefix'")
    end
    
    # Add time estimate
    time_estimate = estimate_completion_time(L, max_weight)
    println("  ⏱️  Estimated time: $time_estimate")
    println("  🔄 Translation symmetries will be identified and deduplicated")
    println()
    
    # Add boundary condition to prefix
    full_prefix = isempty(prefix) ? "periodic" : "$(prefix)_periodic"
    
    # Generate clusters using the same optimized approach as generate_ising_clusters.jl
    println("🚀 Starting single-site cluster enumeration...")
    println("📊 This will involve 5 main steps:")
    println("   1️⃣  Center site loop generation")
    println("   2️⃣  Loop translation and expansion")
    println("   3️⃣  Cluster enumeration for target site")
    println("   4️⃣  Translation-aware deduplication")
    println("   5️⃣  Data saving and analysis")
    println()
    
    start_time = time()
    
    try
        # Step 1: Generate loops supported on center site
        center_site = calculate_center_site(L)
        center_coords = site_to_coords(center_site, L)
        println("📍 Step 1: Generating loops on center site $center_site (coordinates $center_coords)...")
        
        # Use SquareLoopEnumerator for optimized performance
        enumerator_square = SquareLoopEnumerator(L; periodic=true)
        center_loops = find_loops_supported_on_vertex_square(enumerator_square, center_site, max_loop_weight)
        
        println("  Found $(length(center_loops)) loops supported on center site")
        
        # Step 2: Generate all loops via translation
        println("\n🔄 Step 2: Generating all loops via translation...")
        all_loops = generate_all_loops_pbc(center_loops, L)
        
        # Step 3: Build interaction graph and enumerate clusters for target site
        println("\n🔗 Step 3: Building interaction graph and enumerating clusters for site $site...")
        interaction_graph = build_interaction_graph_optimized(all_loops)
        
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
            
            # Create empty data structure
            data = SingleSiteClusterData(
                Cluster[],
                all_loops,
                enumerator_square.base_enumerator.adj_matrix,
                max_weight,
                L,
                site,
                0.0,
                0.0,
                string(now()),
                boundary,
                0,
                0,
                0
            )
            
            save_single_site_cluster_data(data, full_prefix)
            println("\n✅ Empty enumeration completed!")
            return
        end
        
        # DFS cluster enumeration
        site_clusters = dfs_enumerate_clusters_from_supported(all_loops, supported_loop_ids, 
                                                            max_weight, interaction_graph)
        
        # Remove standard redundancies first
        unique_clusters = remove_cluster_redundancies(site_clusters)
        clusters_before_translation = length(unique_clusters)
        
        enum_time = time() - start_time
        
        # Step 4: Apply translation-aware deduplication
        println("\n🔧 Step 4: Translation-aware deduplication...")
        translation_start_time = time()
        
        clusters_after_translation = apply_translation_aware_deduplication(unique_clusters, all_loops, L)
        clusters_after_count = length(clusters_after_translation)
        
        # Count unique canonical forms
        canonical_forms = Set{Tuple}()
        for cluster in clusters_after_translation
            canonical_sig = create_translation_aware_cluster_signature(cluster, all_loops, L)
            push!(canonical_forms, canonical_sig)
        end
        canonical_forms_count = length(canonical_forms)
        
        translation_time = time() - translation_start_time
        
        # Step 5: Create data structure and save
        println("\n💾 Step 5: Creating data structure and saving...")
        total_time = time() - start_time
        
        data = SingleSiteClusterData(
            clusters_after_translation,
            all_loops,
            enumerator_square.base_enumerator.adj_matrix,
            max_weight,
            L,
            site,
            enum_time,
            translation_time,
            string(now()),
            boundary,
            clusters_before_translation,
            clusters_after_count,
            canonical_forms_count
        )
        
        # Save results
        filepath = save_single_site_cluster_data(data, full_prefix)
        
        println("\n✅ Single-site enumeration completed!")
        
        # Success summary
        println("\n🎉 SUCCESS!")
        println("="^40)
        println("✅ Enumeration completed in $(round(total_time, digits=2)) seconds")
        println("✅ Clusters before translation removal: $clusters_before_translation")
        println("✅ Clusters after translation removal: $clusters_after_count")
        println("✅ Unique canonical forms: $canonical_forms_count")
        if clusters_before_translation > 0
            reduction_factor = round(clusters_before_translation / max(1, clusters_after_count), digits=2)
            println("✅ Translation redundancy factor: $(reduction_factor)x")
        end
        println("✅ Data saved successfully!")
        
        # Show file location
        println("📁 Saved to: $(basename(filepath))")
        
        # Quick analysis
        println("\n📊 Quick Analysis:")
        by_weight = Dict{Int, Int}()
        for cluster in clusters_after_translation
            by_weight[cluster.weight] = get(by_weight, cluster.weight, 0) + 1
        end
        
        for weight in sort(collect(keys(by_weight)))
            count = by_weight[weight]
            println("  Weight $weight: $count unique clusters (translation-inequivalent)")
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
    println("  • Note: Only unique clusters under translation symmetry are saved")
end

function show_examples()
    """Show usage examples."""
    println("📚 EXAMPLES:")
    println()
    println("1. Generate clusters for center site on 4×4 periodic lattice, max weight 3:")
    println("   julia generate_ising_clusters_one_site.jl --size 4 --weight 3")
    println()
    println("2. Generate clusters for site 10 on 6×6 lattice with custom prefix:")
    println("   julia generate_ising_clusters_one_site.jl -L 6 -s 10 -w 5 -p my_test")
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