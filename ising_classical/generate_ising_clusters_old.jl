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

using ArgParse
using ProgressMeter

# Helper function to estimate completion time
function estimate_completion_time(L::Int, max_weight::Int)
    """Provide rough time estimates based on lattice size and weight."""
    # Very rough estimates based on complexity
    base_time = L^2 * max_weight^2 * 0.001  # seconds
    
    if L <= 4
        return "< 1 second"
    elseif L <= 6
        return "< 30 seconds"
    elseif L <= 8
        return "1-5 minutes"
    elseif L <= 10
        return "5-30 minutes"
    else
        return "30+ minutes (consider smaller parameters)"
    end
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
    
    # Generate clusters with overall progress tracking
    println("\n🚀 Starting cluster enumeration...")
    println("📊 This will involve 4 main steps:")
    println("   1️⃣  Loop enumeration")
    println("   2️⃣  Interaction graph building")
    println("   3️⃣  Cluster enumeration per site")
    println("   4️⃣  Data saving and analysis")
    println()
    
    start_time = time()
    
    try
        # Track enumeration progress
        println("⏳ Enumeration will begin shortly...")
        println("   (Progress bars will appear for each major step)")
        println()
        
        data = enumerate_all_clusters_all_sites(enumerator, max_weight, max_loop_weight; 
                                               prefix=full_prefix)
        
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