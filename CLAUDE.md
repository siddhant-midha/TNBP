# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a research codebase implementing Tensor Network Belief Propagation (TNBP) for statistical mechanics calculations, primarily focused on the 2D Ising model. The codebase combines belief propagation methods with loop corrections to compute partition functions and free energies.

## Core Architecture

### Main Modules

1. **BP.jl** - Core belief propagation implementation
   - `get_adj_mat()` - Extracts adjacency matrices from tensor networks
   - `get_messages()` - Initializes BP messages between tensor nodes
   - Message passing algorithms for tensor network contraction

2. **Ising2D.jl** - 2D Ising model specific functions
   - `get_ising_tn()` - Generates Ising tensor networks on periodic lattices
   - `periodic_square_lattice()` - Creates periodic square lattice graphs
   - `free_energy()` - Analytical Onsager solution for comparison

3. **boundary_evolution.jl** - PEPS contraction using boundary MPS method
   - `contract_peps_no_phys()` - Contracts 2D PEPS without physical indices
   - Row-by-row contraction using boundary MPS evolution

4. **brute_force.jl** - Loop enumeration and correction calculations
   - Functions for enumerating closed loops of various orders (4th, 6th, 7th, 8th, 9th, 10th)
   - Loop contribution calculations for partition function corrections

### Key Dependencies

The codebase relies on Julia packages loaded in `dependencies.jl`:
- **ITensors.jl** - Core tensor operations and MPS/MPO functionality
- **Graphs.jl** - Graph operations for lattice structures
- **HCubature.jl** - Numerical integration for analytical comparisons
- **Plots.jl** - Visualization of results
- **ProgressMeter.jl** - Progress tracking for long computations

## Working with Jupyter Notebooks

The main analysis is conducted in Jupyter notebooks:

- **loop_correction.ipynb** - Primary notebook implementing loop-corrected BP
- **convergence.ipynb** - Analysis of BP convergence properties
- **2DIsing.ipynb** - 2D Ising model specific calculations
- **operator_dynamics.ipynb** - Time evolution studies

When working with notebooks, use:
```julia
include("dependencies.jl")
include("BP.jl") 
include("Ising2D.jl")
include("boundary_evolution.jl")
include("brute_force.jl")
```

## Common Workflows

### Running Loop-Corrected BP Calculations

Typical workflow for computing free energies with loop corrections:

1. Generate Ising tensor network: `T,indmat = get_ising_tn(g,N,β;h=h)`
2. Initialize and run BP: `messages = message_passing(T,g,messages,indmat)`
3. Compute BP fixed point: `Z_l = get_fixed_point_list(T,g,messages,indmat,N)`
4. Apply loop corrections: `logZCorrection4th()`, `logZCorrection6th()`, etc.
5. Compare with exact Onsager solution and boundary evolution results

### Parameter Scanning

The codebase supports systematic parameter scans:
- β (inverse temperature) scans to study phase transitions
- h (magnetic field) scans for field-dependent behavior
- System size L scaling studies

### Comparison Methods

Three main approaches are implemented for cross-validation:
1. **BP Fixed Point** - Mean field approximation
2. **Loop-Corrected BP** - Systematic corrections up to 10th order
3. **Boundary Evolution** - Exact PEPS contraction using MPS methods
4. **Analytical (Onsager)** - Exact solution for zero field case

## Development Notes

- All tensor operations use ITensors.jl conventions
- Graph structures use Graphs.jl with periodic boundary conditions
- Loop corrections require careful index management between tensor nodes
- Performance-critical sections include message passing and loop enumeration
- Memory usage scales with bond dimensions and system size

## Testing and Validation

Validate implementations by:
- Comparing BP results with exact Onsager solution at h=0
- Checking convergence of loop corrections with increasing order
- Verifying boundary evolution results match analytical predictions
- Testing different system sizes for finite-size effects


given a graph G, a loop is defined as a connected subgraph of G where each vertex has degree >=2. two loops are incompatible if they share at least one vertex. a cluster W is a multiset of loops (so each loop can show up multiple times). given a cluster W with n loops (including multiplicities), the interaction graph G_W contains n vertices corresponding to loops. two vertices are connected if the corresponding loops are incompatible. W is called connected if G_W is a connected graph. the weight of the loop is defined as the number of edges. and the weight of a cluster is defined as the sum of weights of all loops. a cluster is supported on a site if at least one of its loops is supported on that site. 

