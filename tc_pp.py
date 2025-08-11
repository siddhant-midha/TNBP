#!/usr/bin/env python3
import numpy as np
import csv
import os
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict

def load_and_average_data(start_id=1, end_id=100, date="250809", output_filename=None):
    """
    Load CSV data from multiple runs and average over all shots.
    
    Parameters:
    - start_id: Starting ID number (default: 1)
    - end_id: Ending ID number (default: 100)  
    - date: Date folder name (default: "250809")
    - output_filename: Output filename (default: "{start_id}_{end_id}_pp.npz")
    """
    
    if output_filename is None:
        output_filename = f"{start_id}_{end_id}_pp.npz"
    
    data_dir = Path("data") / date
    
    # Dictionary to store data grouped by (L, p) pairs
    grouped_data = defaultdict(list)
    missing_files = []
    successful_files = 0
    
    print(f"Loading CSV files from {data_dir}...")
    
    # Load all CSV files with progress bar
    for file_id in tqdm(range(start_id, end_id + 1), desc="Loading files"):
        csv_file = data_dir / f"{file_id}.csv"
        
        if csv_file.exists():
            try:
                with open(csv_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Convert string values to float
                        L = int(row['L'])
                        p = float(row['p'])
                        key = (L, p)
                        
                        # Store all metrics for this (L, p) combination
                        data_point = {
                            'logical_mean_loops': float(row['logical_mean_loops']),
                            'logical_std_loops': float(row['logical_std_loops']),
                            'failure_mean_loops': float(row['failure_mean_loops']),
                            'failure_std_loops': float(row['failure_std_loops']),
                            'logical_mean_no_loops': float(row['logical_mean_no_loops']),
                            'logical_std_no_loops': float(row['logical_std_no_loops']),
                            'failure_mean_no_loops': float(row['failure_mean_no_loops']),
                            'failure_std_no_loops': float(row['failure_std_no_loops'])
                        }
                        grouped_data[key].append(data_point)
                
                successful_files += 1
                
            except Exception as e:
                print(f"Error loading {csv_file}: {e}")
                missing_files.append(file_id)
        else:
            missing_files.append(file_id)
    
    if not grouped_data:
        raise ValueError(f"No data loaded from {data_dir}")
    
    print(f"Successfully loaded {successful_files} files")
    if missing_files:
        print(f"Missing files for IDs: {missing_files}")
    
    # Convert grouped data to averaged arrays
    print("Computing averages...")
    
    # Sort keys for consistent ordering
    sorted_keys = sorted(grouped_data.keys())
    
    L_values = []
    p_values = []
    metrics_avg = defaultdict(list)
    metrics_std = defaultdict(list)
    
    for key in tqdm(sorted_keys, desc="Computing statistics"):
        L, p = key
        data_points = grouped_data[key]
        
        L_values.append(L)
        p_values.append(p)
        
        # Compute mean and std for each metric across all runs
        for metric in data_points[0].keys():
            values = [point[metric] for point in data_points]
            metrics_avg[metric].append(np.mean(values))
            metrics_std[metric].append(np.std(values, ddof=1) if len(values) > 1 else 0.0)
    
    # Convert lists to numpy arrays
    data_dict = {
        'L': np.array(L_values),
        'p': np.array(p_values),
        'logical_mean_loops_avg': np.array(metrics_avg['logical_mean_loops']),
        'logical_mean_loops_std': np.array(metrics_std['logical_mean_loops']),
        'failure_mean_loops_avg': np.array(metrics_avg['failure_mean_loops']),
        'failure_mean_loops_std': np.array(metrics_std['failure_mean_loops']),
        'logical_mean_no_loops_avg': np.array(metrics_avg['logical_mean_no_loops']),
        'logical_mean_no_loops_std': np.array(metrics_std['logical_mean_no_loops']),
        'failure_mean_no_loops_avg': np.array(metrics_avg['failure_mean_no_loops']),
        'failure_mean_no_loops_std': np.array(metrics_std['failure_mean_no_loops']),
        'num_files_processed': successful_files,
        'missing_files': np.array(missing_files) if missing_files else np.array([]),
        'date': date,
        'id_range': f"{start_id}-{end_id}"
    }
    
    # Save to NPZ file
    print(f"Saving results to {output_filename}...")
    np.savez_compressed(output_filename, **data_dict)
    
    print(f"\nProcessing complete!")
    print(f"Processed {successful_files} files from ID {start_id} to {end_id}")
    print(f"Results saved to: {output_filename}")
    print(f"Data shape: L values: {len(np.unique(L_values))}, p values: {len(np.unique(p_values))}")
    
    return data_dict

if __name__ == "__main__":
    # Default parameters - can be modified as needed
    start_id = 1
    end_id = 100
    
    # Load and process data
    result = load_and_average_data(start_id=start_id, end_id=end_id)
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"Unique L values: {np.unique(result['L'])}")
    print(f"p value range: {result['p'].min():.6f} to {result['p'].max():.6f}")
    print(f"Number of (L,p) combinations: {len(result['L'])}")