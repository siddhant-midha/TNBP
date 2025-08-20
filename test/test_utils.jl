#!/usr/bin/env julia

"""
test_utils.jl

Utility functions for loading test data consistently across all test files.
"""

using Serialization

function load_latest_cluster_file_by_criteria(; size_filter="L11", weight_filter="", boundary_filter="periodic")
    """
    Load the most recent cluster enumeration file matching criteria.
    Returns the data and filename for reference.
    """
    save_dir = "saved_clusters"
    
    if !isdir(save_dir)
        error("No saved_clusters directory found!")
    end
    
    files = readdir(save_dir)
    
    # Filter files based on criteria
    matching_files = filter(files) do f
        # Must be a .jld2 file
        !endswith(f, ".jld2") && return false
        
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
        error("No matching cluster files found for criteria: size=$size_filter, weight=$weight_filter, boundary=$boundary_filter")
    end
    
    # Sort by filename (which includes timestamp) and take the most recent
    latest_file = sort(matching_files)[end]
    filepath = joinpath(save_dir, latest_file)
    
    println("📖 Loading latest cluster data: $(latest_file)")
    
    loaded_data = open(filepath, "r") do io
        deserialize(io)
    end
    
    return loaded_data["data"], latest_file
end

function load_latest_L11_w10_file()
    """Load the most recent L11 w10 periodic file."""
    return load_latest_cluster_file_by_criteria(size_filter="L11", weight_filter="w10", boundary_filter="periodic")
end

function load_latest_L11_file()
    """Load the most recent L11 periodic file (any weight)."""
    return load_latest_cluster_file_by_criteria(size_filter="L11", boundary_filter="periodic")
end

function load_latest_test_file()
    """Load the most recent test file for deduplication testing."""
    save_dir = "saved_clusters"
    files = readdir(save_dir)
    
    # Look for recent test files with fixed deduplication, prioritizing by recency
    test_patterns = ["test_w6_fixed", "test_dedup", "test_w6", "test_"]
    
    for pattern in test_patterns
        matching_files = filter(f -> contains(f, pattern) && 
                                    contains(f, "L11") && 
                                    endswith(f, ".jld2"), files)
        
        if !isempty(matching_files)
            latest_file = sort(matching_files)[end]
            filepath = joinpath(save_dir, latest_file)
            println("📖 Loading latest test data: $(latest_file)")
            
            loaded_data = open(filepath, "r") do io
                deserialize(io)
            end
            return loaded_data["data"], latest_file
        end
    end
    
    error("No suitable test files found for deduplication testing")
end