#!/usr/bin/env python3
import numpy as np
import csv
import os
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict

def load_and_average_data(start_id=1, end_id=100, date="250809", output_filename=None, samples_per_file=1000):
    """
    Load CSV data from multiple runs and compute proper pooled statistics.
    
    Parameters:
    - start_id: Starting ID number (default: 1)
    - end_id: Ending ID number (default: 100)  
    - date: Date folder name (default: "250809")
    - output_filename: Output filename (default: "{start_id}_{end_id}_pp.npz")
    - samples_per_file: Number of samples per CSV file (default: 1000)
    """

    data_dir = Path("data") / date

    if output_filename is None:
        output_filename = f"{start_id}_{end_id}_pp.npz"
    output_filename = data_dir / output_filename
    
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
    metrics_pooled_std = defaultdict(list)  # Pooled standard error
    num_files_per_key = defaultdict(int)
    
    for key in tqdm(sorted_keys, desc="Computing statistics"):
        L, p = key
        data_points = grouped_data[key]
        
        L_values.append(L)
        p_values.append(p)
        num_files_per_key[key] = len(data_points)
        
        # Compute pooled mean and standard error for each metric
        for metric_base in ['logical_mean_loops', 'failure_mean_loops', 'logical_mean_no_loops', 'failure_mean_no_loops']:
            mean_key = metric_base
            std_key = metric_base.replace('_mean_', '_std_')
            
            # Extract means and standard deviations from all files
            means = [point[mean_key] for point in data_points]
            stds = [point[std_key] for point in data_points]
            
            # Compute pooled mean (weighted average)
            pooled_mean = np.mean(means)
            
            # Compute pooled variance using the formula for combining samples
            # Var_pooled = (1/N_total) * Σ[(n_i - 1) * s_i^2 + n_i * (x̄_i - x̄_pooled)^2]
            total_samples = len(data_points) * samples_per_file
            pooled_variance = 0.0
            
            for i, (mean_i, std_i) in enumerate(zip(means, stds)):
                # Variance from within-group variation
                within_var = (samples_per_file - 1) * (std_i ** 2)
                # Variance from between-group variation  
                between_var = samples_per_file * ((mean_i - pooled_mean) ** 2)
                pooled_variance += (within_var + between_var)
            
            pooled_variance /= (total_samples - 1)  # Bessel's correction
            pooled_std = np.sqrt(pooled_variance)
            
            # Standard error of the pooled mean
            standard_error = pooled_std / np.sqrt(total_samples)
            
            metrics_avg[mean_key].append(pooled_mean)
            metrics_pooled_std[mean_key].append(standard_error)
    
    # Convert lists to numpy arrays
    data_dict = {
        'L': np.array(L_values),
        'p': np.array(p_values),
        'logical_mean_loops_avg': np.array(metrics_avg['logical_mean_loops']),
        'logical_mean_loops_std': np.array(metrics_pooled_std['logical_mean_loops']),
        'failure_mean_loops_avg': np.array(metrics_avg['failure_mean_loops']),
        'failure_mean_loops_std': np.array(metrics_pooled_std['failure_mean_loops']),
        'logical_mean_no_loops_avg': np.array(metrics_avg['logical_mean_no_loops']),
        'logical_mean_no_loops_std': np.array(metrics_pooled_std['logical_mean_no_loops']),
        'failure_mean_no_loops_avg': np.array(metrics_avg['failure_mean_no_loops']),
        'failure_mean_no_loops_std': np.array(metrics_pooled_std['failure_mean_no_loops']),
        'num_files_processed': successful_files,
        'missing_files': np.array(missing_files) if missing_files else np.array([]),
        'date': date,
        'id_range': f"{start_id}-{end_id}",
        'total_samples_per_point': len(grouped_data[sorted_keys[0]]) * samples_per_file if sorted_keys else 0
    }
    
    # Save to NPZ file
    print(f"Saving results to {output_filename}...")
    np.savez_compressed(output_filename, **data_dict)
    
    print(f"\nProcessing complete!")
    print(f"Processed {successful_files} files from ID {start_id} to {end_id}")
    print(f"Results saved to: {output_filename}")
    print(f"Data shape: L values: {len(np.unique(L_values))}, p values: {len(np.unique(p_values))}")
    print(f"Total samples per (L,p) point: {len(grouped_data[sorted_keys[0]]) * samples_per_file if sorted_keys else 0}")
    
    return data_dict

if __name__ == "__main__":
    # Default parameters - can be modified as needed
    start_id = 201
    end_id = 600
    samples_per_file = 1000  # Adjust this to match your actual number of samples per CSV file
    
    # Load and process data
    result = load_and_average_data(start_id=start_id, end_id=end_id, samples_per_file=samples_per_file)
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"Unique L values: {np.unique(result['L'])}")
    print(f"p value range: {result['p'].min():.6f} to {result['p'].max():.6f}")
    print(f"Number of (L,p) combinations: {len(result['L'])}")
    print(f"Total samples per data point: {result['total_samples_per_point']}")
    print(f"Error bars now represent standard error of {result['total_samples_per_point']} pooled samples")