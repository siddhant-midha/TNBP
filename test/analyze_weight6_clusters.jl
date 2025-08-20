#!/usr/bin/env julia

"""
analyze_weight6_clusters.jl

Analyze the structure of weight-6 clusters to understand why there are 442 instead of expected 242.
"""

include("../dependencies.jl")
include("../functions/ClusterEnumeration.jl")

using Serialization

function load_latest_w6_file()
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    matching_files = filter(f -> contains(f, "L11") && contains(f, "w6") && contains(f, "test_w6") && endswith(f, ".jld2"), files)
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    println("📖 Loading: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    return loaded_data["data"]
end

function analyze_weight6_cluster_types()
    println("🔍 Analyzing Weight-6 Cluster Types")
    println("="^50)
    
    data = load_latest_w6_file()
    L = data.lattice_size
    
    # Collect all weight-6 clusters
    weight6_clusters = []
    for (site, clusters) in data.clusters_by_site
        for cluster in clusters
            if cluster.weight == 6
                push!(weight6_clusters, (site, cluster))
            end
        end
    end
    
    println("📊 Found $(length(weight6_clusters)) weight-6 clusters total")
    println("📊 Expected: 121*2 = 242 clusters")
    
    # Analyze cluster structures
    structure_counts = Dict()
    
    for (site, cluster) in weight6_clusters
        # Analyze the structure
        n_loops = length(cluster.loop_ids)
        multiplicities = sort(collect(values(cluster.multiplicities)))
        
        # Get loop weights
        loop_weights = []
        for loop_id in cluster.loop_ids
            loop = data.all_loops[loop_id]
            push!(loop_weights, loop.weight)
        end
        sort!(loop_weights)
        
        structure_key = (n_loops, multiplicities, loop_weights)
        structure_counts[structure_key] = get(structure_counts, structure_key, 0) + 1
    end
    
    println("\n📋 Weight-6 cluster structures:")
    total_accounted = 0
    for (structure, count) in sort(collect(structure_counts), by=x->x[2], rev=true)
        n_loops, multiplicities, loop_weights = structure
        println("  $count clusters: $n_loops loops, multiplicities=$multiplicities, weights=$loop_weights")
        
        # Verify weight calculation
        total_weight = sum(loop_weights[i] * multiplicities[i] for i in 1:length(loop_weights))
        if total_weight != 6
            println("    ❌ ERROR: Total weight = $total_weight ≠ 6")
        else
            total_accounted += count
        end
    end
    
    println("\n📊 Total accounted: $total_accounted / $(length(weight6_clusters))")
    
    # Look at specific examples of each type
    println("\n🔍 Detailed examples by structure type:")
    
    structure_examples = Dict()
    for (site, cluster) in weight6_clusters
        n_loops = length(cluster.loop_ids)
        multiplicities = sort(collect(values(cluster.multiplicities)))
        loop_weights = [data.all_loops[loop_id].weight for loop_id in cluster.loop_ids]
        sort!(loop_weights)
        
        structure_key = (n_loops, multiplicities, loop_weights)
        
        if !haskey(structure_examples, structure_key)
            structure_examples[structure_key] = []
        end
        if length(structure_examples[structure_key]) < 3  # Keep first 3 examples
            push!(structure_examples[structure_key], (site, cluster))
        end
    end
    
    for (structure, examples) in structure_examples
        n_loops, multiplicities, loop_weights = structure
        count = structure_counts[structure]
        println("\n  Structure: $n_loops loops, mults=$multiplicities, weights=$loop_weights ($count total)")
        
        for (i, (site, cluster)) in enumerate(examples)
            println("    Example $i (site $site):")
            println("      Loop IDs: $(cluster.loop_ids)")
            println("      Multiplicities: $(cluster.multiplicities)")
            
            # Show loop details
            for loop_id in cluster.loop_ids
                loop = data.all_loops[loop_id]
                mult = cluster.multiplicities[loop_id]
                println("        Loop $loop_id (×$mult): weight=$(loop.weight), vertices=$(loop.vertices)")
            end
        end
    end
    
    # Analyze the loop types that exist
    println("\n🔍 Loop type analysis:")
    loop_types = Dict()
    
    for loop in data.all_loops
        loop_types[loop.weight] = get(loop_types, loop.weight, 0) + 1
    end
    
    println("  Loop distribution:")
    for weight in sort(collect(keys(loop_types)))
        count = loop_types[weight]
        println("    Weight $weight: $count loops")
    end
    
    # Check if there are weight-6 loops (single loops with weight 6)
    weight6_loops = [loop for loop in data.all_loops if loop.weight == 6]
    println("\n  Weight-6 loops: $(length(weight6_loops))")
    
    if length(weight6_loops) > 0
        println("  Expected weight-6 single-loop clusters: $(length(weight6_loops))")
        
        # Count actual single weight-6 loop clusters
        single_w6_clusters = 0
        for (site, cluster) in weight6_clusters
            if length(cluster.loop_ids) == 1 && 
               data.all_loops[cluster.loop_ids[1]].weight == 6 &&
               cluster.multiplicities[cluster.loop_ids[1]] == 1
                single_w6_clusters += 1
            end
        end
        println("  Actual single weight-6 loop clusters: $single_w6_clusters")
    end
    
    # Check for clusters made from multiple loops
    println("\n🔍 Multi-loop cluster analysis:")
    multi_loop_types = Dict()
    
    for (site, cluster) in weight6_clusters
        if length(cluster.loop_ids) > 1
            # Get composition
            composition = []
            for loop_id in sort(cluster.loop_ids)
                loop = data.all_loops[loop_id]
                mult = cluster.multiplicities[loop_id]
                push!(composition, (loop.weight, mult))
            end
            sort!(composition)
            
            comp_key = tuple(composition...)
            multi_loop_types[comp_key] = get(multi_loop_types, comp_key, 0) + 1
        end
    end
    
    println("  Multi-loop cluster compositions:")
    for (composition, count) in sort(collect(multi_loop_types), by=x->x[2], rev=true)
        total_weight = sum(weight * mult for (weight, mult) in composition)
        println("    $count clusters: $composition (total weight: $total_weight)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_weight6_cluster_types()
end