# Code Structure Overview

This directory contains complete implementations for enumerating loops and clusters on periodic square lattices in both Julia and Python. Here's what each file contains:

## 🎯 **Current Recommended Implementation: Julia**

### 1. `LoopEnumeration.jl` ⭐
**Purpose:** Main Julia loop enumeration (degree ≥ 2, allows intersections)
**Main Functions:**
- `LoopEnumerator` struct - handles complex loop structures
- `find_loops_supported_on_vertex()` - BFS-based expansion algorithm
- `is_valid_loop()` - validates degree ≥ 2 constraint
- `create_periodic_square_lattice()` - generates lattice adjacency matrix

### 2. `ClusterEnumeration.jl` 
**Purpose:** Julia cluster enumeration (multisets of loops)
**Main Functions:**
- `ClusterEnumerator` struct - finds connected clusters of loops
- `enumerate_connected_clusters()` - main cluster enumeration function
- `build_interaction_graph()` - creates incompatibility graph

### 3. `test_loop_enumeration.jl` ⭐
**Purpose:** Comprehensive Julia test suite for loop enumeration
**Test Functions:**
- `test_10x10_loop_enumeration()` - main 10×10 lattice test
- `test_smaller_lattices()` - pattern analysis across sizes
- `verify_loop_connectivity()` - validates found loops

### 4. `test_cluster_enumeration.jl`
**Purpose:** Julia cluster enumeration testing
**Test Functions:**
- `test_10x10_cluster_enumeration()` - complete cluster test
- `analyze_cluster_results()` - detailed result analysis

### 5. `simple_test.jl` ⭐
**Purpose:** Quick demonstration of working Julia implementation
**Shows:** 28 loops found correctly on 10×10 lattice in ~0.4 seconds

## 📁 **Python Implementation (Legacy - in ./python/ folder)**

### Core Python Files (in ./python/):

### 1. `python/optimized_loop_enumeration.py` ⭐
**Purpose:** Python loop enumeration (degree ≥ 2, allows intersections)
**Main Functions:**
- `OptimizedLoopEnumerator` class - handles complex loop structures
- `find_loops_supported_on_vertex()` - BFS-based expansion algorithm

### 2. `python/cluster_enumeration_complete.py` ⭐
**Purpose:** Complete Python cluster enumeration (multisets of loops)
**Main Functions:**
- `ClusterEnumerator` class - finds connected clusters of loops
- `enumerate_connected_clusters()` - main cluster enumeration function

### 3. `python/test_10x10_corrected.py` ⭐
**Purpose:** Tests for corrected Python loop enumeration
**Test Functions:**
- `test_10x10_corrected()` - tests degree ≥ 2 loops on 10×10 lattice

### 4. `python/test_cluster_10x10.py` ⭐
**Purpose:** Main Python cluster enumeration testing
**Test Functions:**
- `test_10x10_cluster_enumeration()` - complete cluster test

### Other Python Files (Development History):
- `python/loop_enumeration.py` - Original simple cycle enumeration
- `python/fixed_loop_enumeration.py` - Improved simple cycles  
- `python/check_redundancy.py` - Duplicate verification
- `python/debug_*.py` - Various debugging utilities
- `python/analyze_expected_loops.py` - Theoretical analysis
- `python/manual_verification.py` - Manual counting verification

## 📋 Documentation Files

### 14. `cluster_enumeration.jl`
**Purpose:** Documentation file explaining the Python implementation
**Content:** Usage instructions and algorithm description

### 15. `CLAUDE.md`
**Purpose:** Development guide for the entire TNBP project
**Content:** Project overview, architecture, workflows

### 16. `CODE_OVERVIEW.md` (this file)
**Purpose:** Explains the structure of all code files

## 🎯 Which Files to Use

### For Loop Enumeration:
- **Use:** `optimized_loop_enumeration.py` 
- **Why:** Correctly handles degree ≥ 2 constraint, finds all 28 loops

### For Cluster Enumeration:  
- **Use:** `cluster_enumeration_complete.py`
- **Why:** Complete implementation of multiset cluster enumeration

### For Testing:
- **Quick tests:** `test_10x10_corrected.py`
- **Full cluster tests:** `test_cluster_10x10.py`
- **Verification:** `check_redundancy.py`

### For Understanding:
- **Algorithm development:** Read files in order: `loop_enumeration.py` → `fixed_loop_enumeration.py` → `optimized_loop_enumeration.py`
- **Problem evolution:** See how the definition of "loop" evolved from simple cycles to general degree ≥ 2 structures

## ⚡ Quick Start

### Julia (Recommended - Fast):
```bash
# Quick test (finds 28 loops in ~0.4s)
julia simple_test.jl

# Full loop enumeration test
julia test_loop_enumeration.jl

# Cluster enumeration test  
julia test_cluster_enumeration.jl
```

### Python (Legacy - in ./python/ folder):
```bash
# Test loop enumeration (finds 28 loops in ~3s)
python python/test_10x10_corrected.py

# Test cluster enumeration (finds 28 clusters)  
python python/test_cluster_10x10.py

# Verify no duplicates
python python/check_redundancy.py
```

## 🔄 Algorithm Evolution

1. **Simple Cycles** (`loop_enumeration.py`) → 16 loops (degree = 2 only)
2. **Fixed Simple Cycles** (`fixed_loop_enumeration.py`) → 16 loops (better algorithm)  
3. **General Loops** (`optimized_loop_enumeration.py`) → 28 loops (degree ≥ 2) ✅
4. **Clusters** (`cluster_enumeration_complete.py`) → 28 clusters (multisets) ✅

The evolution shows how understanding the correct problem definition (degree ≥ 2) was crucial for getting the expected result of 28.