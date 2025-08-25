#!/usr/bin/env julia

"""
generate_ising_loops_one_site.jl

Generate and save loop enumerations for a single site on square lattice Ising models.
Supports both periodic and open boundary conditions.
Includes translational symmetry removal - only keeps one representative from each equivalence class.

Usage:
    julia generate_ising_loops_one_site.jl --size 6 --weight 8 --site 1 --boundary periodic --prefix test
    julia generate_ising_loops_one_site.jl --help
"""

include("../dependencies.jl")
include("../functions/LoopEnumeration.jl")
include("../functions/loop_enumeration_square.jl")

using ArgParse
using ProgressMeter
using Dates
using Serialization

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
function estimate_completion_time_loops(L::Int, max_weight::Int)
    """Provide rough time estimates based on lattice size and weight for loop enumeration."""
    if L <= 6
        return "< 1 second"
    elseif L <= 8
        return "< 10 seconds"
    elseif L <= 10
        return "10-60 seconds"
    elseif L <= 12
        return "1-5 minutes"
    else
        return "5+ minutes (consider smaller parameters)"
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

# Translation functions for periodic lattices
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

function canonical_loop_under_translation(loop::Loop, L::Int)
    """Find the canonical representative of a loop under translation symmetry."""
    canonical_loop = loop
    canonical_signature = canonical_loop_representation(canonical_loop)
    
    # Try all possible translations
    for di in 0:(L-1)
        for dj in 0:(L-1)
            if di == 0 && dj == 0
                continue  # Skip identity translation
            end
            
            translated_loop = translate_loop(loop, di, dj, L)
            translated_signature = canonical_loop_representation(translated_loop)
            
            # Keep the lexicographically smallest signature
            if translated_signature < canonical_signature
                canonical_loop = translated_loop
                canonical_signature = translated_signature
            end
        end
    end
    
    return canonical_loop
end

function remove_translational_loop_redundancies(loops::Vector{Loop}, L::Int)
    """Remove loops that are related by translation, keeping only canonical representatives."""
    println("🔄 Removing translational redundancies from loops...")
    flush(stdout)
    
    unique_loops = Loop[]
    seen_signatures = Set{Tuple{Vector{Int}, Vector{Tuple{Int,Int}}, Int}}()
    
    progress = Progress(length(loops), dt=0.1, desc="Removing redundancies: ", color=:orange, barlen=50)
    
    for loop in loops
        # Find the canonical representative under translation
        canonical_loop = canonical_loop_under_translation(loop, L)
        canonical_signature = canonical_loop_representation(canonical_loop)
        
        if !(canonical_signature in seen_signatures)
            push!(seen_signatures, canonical_signature)
            push!(unique_loops, canonical_loop)
        end
        
        next!(progress, showvalues = [("Unique", "$(length(unique_loops))"), ("Total", "$(length(loops))")])
    end
    
    finish!(progress)
    
    println("  Original loops: $(length(loops))")
    println("  Unique loops (after translation): $(length(unique_loops))")
    if length(loops) > 0
        println("  Redundancy factor: $(round(length(loops) / length(unique_loops), digits=2))x")
    end
    
    return unique_loops
end

# Data structure for single-site loop enumeration results
struct SingleSiteLoopData
    """Data structure to save single-site loop enumeration results."""
    loops::Vector{Loop}
    adj_matrix::Matrix{Int}
    max_weight::Int
    lattice_size::Int
    site::Int
    enumeration_time::Float64
    translation_removal_time::Float64
    timestamp::String
    boundary_condition::String
    loops_before_translation_removal::Int
    loops_after_translation_removal::Int
end

function enumerate_loops_one_site_with_translation_removal(site::Int, max_weight::Int, L::Int, boundary::String)
    """Enumerate loops for one site and remove translational redundancies."""
    
    println("🚀 Starting single-site loop enumeration...")
    println("📋 Target site: $site")
    println("📋 Max loop weight: $max_weight")
    println("📋 Boundary condition: $boundary")
    println()
    
    start_time = time()
    
    # Step 1: Enumerate loops using original LoopEnumerator (working implementation)
    println("📊 Step 1: Enumerating loops supported on site $site...")
    flush(stdout)
    
    # Create adjacency matrix and enumerator
    println("  🏗️  Setting up lattice and enumerator...")
    setup_progress = Progress(3, dt=0.5, desc="Setup: ", color=:blue, barlen=30)
    
    next!(setup_progress, showvalues = [("Step", "Initializing square enumerator")])
    enumerator = SquareLoopEnumerator(L; periodic=(boundary == "periodic"))
    adj_matrix = enumerator.base_enumerator.adj_matrix  # Get adjacency matrix for data structure
    
    next!(setup_progress, showvalues = [("Step", "Setup complete")])
    finish!(setup_progress)
    
    # Find loops supported on target site with real-time progress indicator
    println("  🔍 Enumerating loops...")
    println("    📍 Target site: $site ($(site_to_coords(site, L)))")
    println("    ⚖️  Max weight: $max_weight")
    
    # Create a progress bar that updates during enumeration
    enum_progress = Progress(1, dt=0.1, desc="Loop enumeration: ", color=:green, barlen=40)
    enum_start = time()
    
    # Run enumeration in a separate task
    result_channel = Channel{Vector{Loop}}(1)
    
    # Start the enumeration in a separate task
    enum_task = @async begin
        loops = find_loops_supported_on_vertex_square(enumerator, site, max_weight)
        put!(result_channel, loops)
    end
    
    # Animate the progress bar while waiting
    site_loops = nothing
    counter = 0
    while !isready(result_channel)
        sleep(0.2)  # Check every 200ms
        counter += 1
        elapsed = round(time() - enum_start, digits=1)
        
        # Show a spinning progress indicator
        spinner = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
        spinner_char = spinner[mod(counter, length(spinner)) + 1]
        
        # Update progress bar with current status
        next!(enum_progress, showvalues = [
            ("Status", "$spinner_char Searching loops"), 
            ("Elapsed", "$(elapsed)s")
        ])
    end
    
    # Get the result
    site_loops = take!(result_channel)
    wait(enum_task)
    
    enum_duration = time() - enum_start
    
    # Complete the progress bar with results
    finish!(enum_progress, showvalues = [
        ("Status", "✅ Complete"), 
        ("Found", "$(length(site_loops)) loops"), 
        ("Time", "$(round(enum_duration, digits=2))s")
    ])
    println()
    loops_before_translation = length(site_loops)
    
    println("📊 Found $(loops_before_translation) loops supported on site $site")
    
    enum_time = time() - start_time
    
    # Step 2: Remove translational redundancies (only for periodic boundary conditions)
    println("\n🔄 Step 2: Removing translational redundancies...")
    translation_start_time = time()
    
    if boundary == "periodic" && !isempty(site_loops)
        loops_after_translation = remove_translational_loop_redundancies(site_loops, L)
    else
        loops_after_translation = site_loops
        if boundary == "open"
            println("  Open boundaries: translation removal not applicable")
        end
    end
    
    loops_after_count = length(loops_after_translation)
    translation_time = time() - translation_start_time
    
    return loops_before_translation, loops_after_translation, loops_after_count, 
           enum_time, translation_time, adj_matrix
end

function save_single_site_loop_data(data::SingleSiteLoopData, prefix::String = "")
    """Save single-site loop enumeration data to file."""
    
    # Create directory if it doesn't exist
    save_dir = "saved_loops"
    if !isdir(save_dir)
        mkpath(save_dir)
        println("📁 Created directory: $save_dir")
    end
    
    # Generate filename
    timestamp = replace(data.timestamp, ":" => "-", "." => "-")
    base_name = "single_site_loops_L$(data.lattice_size)_site$(data.site)_w$(data.max_weight)_$(data.boundary_condition)_$(timestamp)"
    if !isempty(prefix)
        base_name = "$(prefix)_$(base_name)"
    end
    
    filepath = joinpath(save_dir, "$(base_name).jld2")
    
    # Save using Julia's serialization
    println("💾 Saving single-site loop data to: $(filepath)")
    
    # Show progress for save operation
    save_progress = Progress(2, dt=0.2, desc="Saving: ", color=:cyan, barlen=30)
    
    next!(save_progress, showvalues = [("Step", "Preparing data")])
    
    # Create summary for quick inspection
    summary = Dict(
        "lattice_size" => data.lattice_size,
        "site" => data.site,
        "max_weight" => data.max_weight,
        "boundary_condition" => data.boundary_condition,
        "loops_before_translation_removal" => data.loops_before_translation_removal,
        "loops_after_translation_removal" => data.loops_after_translation_removal,
        "redundancy_factor" => data.loops_before_translation_removal / max(1, data.loops_after_translation_removal),
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
    
    next!(save_progress, showvalues = [("Step", "Writing to disk")])
    
    open(filepath, "w") do io
        serialize(io, save_data)
    end
    
    finish!(save_progress)
    println("✅ Saved successfully!")
    println("   Summary: $summary")
    
    return filepath
end

function analyze_loop_statistics(loops::Vector{Loop})
    """Analyze and display statistics about the enumerated loops."""
    println("\n📊 Loop Statistics:")
    
    if isempty(loops)
        println("  No loops found.")
        return
    end
    
    # Group by weight
    by_weight = Dict{Int, Vector{Loop}}()
    for loop in loops
        if !haskey(by_weight, loop.weight)
            by_weight[loop.weight] = Loop[]
        end
        push!(by_weight[loop.weight], loop)
    end
    
    println("Weight | Count | Examples")
    println("-------|-------|----------")
    
    for weight in sort(collect(keys(by_weight)))
        loops_of_weight = by_weight[weight]
        count = length(loops_of_weight)
        
        # Show examples
        examples = String[]
        for i in 1:min(2, count)
            loop = loops_of_weight[i]
            vertices = loop.vertices
            
            # Count degree of each vertex in this loop
            degree = Dict{Int, Int}()
            for v in vertices
                degree[v] = 0
            end
            for (u, v) in loop.edges
                degree[u] += 1
                degree[v] += 1
            end
            
            # Find vertices with degree > 2
            high_degree = [v for v in vertices if degree[v] > 2]
            
            if !isempty(high_degree)
                push!(examples, "V:$(length(vertices)), deg>2:$(length(high_degree))")
            else
                push!(examples, "V:$(length(vertices)), all deg=2")
            end
        end
        
        examples_str = join(examples, "; ")
        if count > 2
            examples_str *= "; +$(count-2) more"
        end
        
        println("  $(lpad(weight, 2))   | $(lpad(count, 4))  | $examples_str")
    end
    
    # Additional statistics
    total_vertices = sum(length(loop.vertices) for loop in loops)
    avg_vertices = round(total_vertices / length(loops), digits=2)
    
    println("\nAdditional Statistics:")
    println("  Total loops: $(length(loops))")
    println("  Average vertices per loop: $avg_vertices")
    println("  Weight range: $(minimum(loop.weight for loop in loops)) - $(maximum(loop.weight for loop in loops))")
end

function list_saved_loop_files()
    """List all saved loop files with basic information."""
    save_dir = "saved_loops"
    
    if !isdir(save_dir)
        println("📂 No saved_loops directory found.")
        return
    end
    
    files = filter(x -> endswith(x, ".jld2"), readdir(save_dir))
    
    if isempty(files)
        println("📂 No saved loop files found in $save_dir/")
        return
    end
    
    println("📂 Saved Loop Files in $save_dir/:")
    println("="^60)
    
    for file in sort(files)
        filepath = joinpath(save_dir, file)
        try
            loaded_data = open(filepath, "r") do io
                deserialize(io)
            end
            
            if haskey(loaded_data, "summary")
                summary = loaded_data["summary"]
                L = summary["lattice_size"]
                site = summary["site"]
                weight = summary["max_weight"]
                boundary = summary["boundary_condition"]
                loops = summary["loops_after_translation_removal"]
                time_str = round(summary["total_time_seconds"], digits=3)
                
                println("📄 $file")
                println("   $(L)×$(L) lattice, site $site, max weight $weight, $boundary BC")
                println("   $(loops) unique loops, $(time_str)s enumeration time")
                println()
            else
                println("📄 $file (invalid format)")
            end
        catch e
            println("📄 $file (error reading: $e)")
        end
    end
end

function parse_commandline()
    s = ArgParseSettings(
        description = "Generate loops for single site on square lattice Ising model",
        epilog = "Example: julia generate_ising_loops_one_site.jl --size 6 --weight 8 --site 1 --boundary periodic"
    )

    @add_arg_table! s begin
        "--size", "-L"
            help = "Linear lattice size (L×L lattice)"
            arg_type = Int
            default = 6
            
        "--weight", "-w"
            help = "Maximum loop weight to enumerate"
            arg_type = Int
            default = 8
            
        "--site", "-s"
            help = "Target site (1-indexed)"
            arg_type = Int
            default = 1
            
        "--boundary", "-b"
            help = "Boundary conditions: 'periodic' or 'open'"
            arg_type = String
            default = "periodic"
            
        "--prefix", "-p"
            help = "Optional prefix for saved files"
            arg_type = String
            default = ""
            
        "--list"
            help = "List existing saved loop files and exit"
            action = :store_true
            
        "--analyze"
            help = "Analyze a saved loop file"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    
    println("🔗 Single-Site Ising Loop Generator")
    println("="^50)
    
    # Handle special actions
    if args["list"]
        list_saved_loop_files()
        return
    end
    
    if !isempty(args["analyze"])
        filepath = args["analyze"]
        if !isfile(filepath)
            # Try looking in saved_loops directory
            filepath = joinpath("saved_loops", filepath)
            if !isfile(filepath)
                println("❌ File not found: $(args["analyze"])")
                return
            end
        end
        
        println("📊 Analyzing: $filepath")
        loaded_data = open(filepath, "r") do io
            deserialize(io)
        end
        
        if haskey(loaded_data, "summary") && haskey(loaded_data, "data")
            summary = loaded_data["summary"]
            data = loaded_data["data"]
            println("✅ Loaded successfully!")
            println("   Summary: $summary")
            
            # Analyze the loops
            analyze_loop_statistics(data.loops)
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
    println("  Max loop weight: $max_weight")
    if !isempty(prefix)
        println("  File prefix: '$prefix'")
    end
    
    # Add time estimate
    time_estimate = estimate_completion_time_loops(L, max_weight)
    println("  ⏱️  Estimated time: $time_estimate")
    
    if boundary == "periodic"
        println("  🔄 Translation symmetries will be removed")
    else
        println("  ⚠️  Open boundaries: translation removal not applicable")
    end
    println()
    
    # Note: Lattice creation is now handled by SquareLoopEnumerator
    println("🏗️  Square lattice will be created during enumeration with $boundary boundary conditions...")
    
    # Add boundary condition to prefix
    full_prefix = isempty(prefix) ? boundary : "$(prefix)_$(boundary)"
    
    start_time = time()
    
    # Overall progress tracking
    total_steps = boundary == "periodic" ? 4 : 3  # Setup, enumerate, translate (if periodic), save
    overall_progress = Progress(total_steps, dt=1.0, desc="Overall progress: ", color=:cyan, barlen=50)
    
    try
        # Step 1: Setup and enumeration
        next!(overall_progress, showvalues = [("Stage", "Loop enumeration")])
        
        # Enumerate loops for the single site
        loops_before_count, loops_after_translation, loops_after_count, 
        enum_time, translation_time, adj_matrix = enumerate_loops_one_site_with_translation_removal(
            site, max_weight, L, boundary)
        
        # Step 2: Translation removal (if applicable)
        if boundary == "periodic"
            next!(overall_progress, showvalues = [("Stage", "Translation removal complete")])
        end
        
        # Step 3: Creating data structure
        next!(overall_progress, showvalues = [("Stage", "Creating data structure")])
        
        total_time = time() - start_time
        
        # Create data structure
        data = SingleSiteLoopData(
            loops_after_translation,
            adj_matrix,
            max_weight,
            L,
            site,
            enum_time,
            translation_time,
            string(now()),
            boundary,
            loops_before_count,
            loops_after_count
        )
        
        # Step 4: Save results
        next!(overall_progress, showvalues = [("Stage", "Saving results")])
        filepath = save_single_site_loop_data(data, full_prefix)
        
        # Complete overall progress
        finish!(overall_progress)
        
        # Success summary
        println("\n🎉 SUCCESS!")
        println("="^40)
        println("✅ Single-site loop enumeration completed in $(round(total_time, digits=2)) seconds")
        println("✅ Loops before translation removal: $loops_before_count")
        println("✅ Loops after translation removal: $loops_after_count")
        if loops_before_count > 0
            reduction_factor = round(loops_before_count / max(1, loops_after_count), digits=2)
            println("✅ Redundancy reduction factor: $(reduction_factor)x")
        end
        println("✅ Data saved successfully!")
        println("📁 Saved to: $(basename(filepath))")
        
        # Loop analysis
        analyze_loop_statistics(loops_after_translation)
        
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
    println("  • Compare with cluster enumeration using generate_ising_clusters_one_site.jl")
end

function show_examples()
    """Show usage examples."""
    println("📚 EXAMPLES:")
    println()
    println("1. Generate loops for site 1 on 4×4 periodic lattice, max weight 6:")
    println("   julia generate_ising_loops_one_site.jl --size 4 --site 1 --weight 6 --boundary periodic")
    println()
    println("2. Generate loops for site 10 on 6×6 open lattice with custom prefix:")
    println("   julia generate_ising_loops_one_site.jl -L 6 -s 10 -w 8 -b open -p my_test")
    println()
    println("3. List all saved loop files:")
    println("   julia generate_ising_loops_one_site.jl --list")
    println()
    println("4. Analyze a saved file:")
    println("   julia generate_ising_loops_one_site.jl --analyze saved_loops/single_site_*.jld2")
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