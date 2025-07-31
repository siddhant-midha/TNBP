# Function Reference Guide

## 🎯 Main Classes and Their Key Methods

### `OptimizedLoopEnumerator` (optimized_loop_enumeration.py)
The core class for finding loops with degree ≥ 2.

```python
enumerator = OptimizedLoopEnumerator(adjacency_matrix)
loops = enumerator.find_loops_supported_on_vertex(vertex, max_weight)
```

**Key Methods:**
- `find_loops_supported_on_vertex(vertex, max_weight)` → List[Dict]
  - Finds all loops containing `vertex` with weight ≤ `max_weight`
  - Returns list of dicts with keys: 'vertices', 'edges', 'weight'

- `_find_loops_bfs(start_vertex, max_weight, found_loops)` 
  - BFS-based loop discovery starting from a vertex
  - Grows subgraphs incrementally

- `_is_valid_loop(vertices, edges)` → bool
  - Validates that all vertices have degree ≥ 2
  - Checks connectivity of the subgraph

### `ClusterEnumerator` (cluster_enumeration_complete.py)  
The main class for finding connected clusters of loops.

```python
enumerator = ClusterEnumerator(adjacency_matrix)
clusters = enumerator.enumerate_connected_clusters(site, max_weight, max_loop_weight)
```

**Key Methods:**
- `enumerate_connected_clusters(site, max_weight, max_loop_weight)` → List[Dict]
  - Main function - finds all connected clusters supported on `site`
  - Returns clusters with keys: 'loop_ids', 'multiplicities', 'weight'

- `_find_all_loops_in_graph(max_weight)` → List[Dict]
  - Discovers all unique loops in the entire graph
  - Removes duplicates using canonical representations

- `_build_interaction_graph(loops)` → Dict[int, List[int]]
  - Creates adjacency list where loops are connected if incompatible
  - Two loops are incompatible if they share vertices

- `_enumerate_clusters_up_to_weight(...)` → List[Dict]
  - Generates all possible multisets of loops up to total weight
  - Uses dynamic programming for efficient enumeration

- `_filter_connected_clusters(clusters, interaction_graph)` → List[Dict]
  - Keeps only clusters whose interaction graph is connected
  - Uses BFS to check connectivity

## 🛠️ Utility Functions

### Lattice Creation
```python
from optimized_loop_enumeration import create_periodic_square_lattice
adj_matrix = create_periodic_square_lattice(L)  # Creates L×L torus
```

### Redundancy Checking
```python
from check_redundancy import remove_duplicates_improved
unique_loops = remove_duplicates_improved(loops)
```

## 🧪 Test Functions

### Basic Loop Testing
```python
# From test_10x10_corrected.py
def test_10x10_corrected():
    """Tests loop enumeration on 10×10 lattice, expects 28 loops"""
    
# From check_redundancy.py  
def analyze_potential_duplicates():
    """Comprehensive duplicate analysis with detailed reporting"""
```

### Cluster Testing
```python
# From test_cluster_10x10.py
def test_10x10_cluster_enumeration():
    """Tests cluster enumeration on 10×10 lattice, expects 28 clusters"""
    
def analyze_cluster_results(clusters, max_weight):
    """Detailed analysis of cluster composition and structure"""
```

## 📊 Data Structures

### Loop Representation
```python
loop = {
    'vertices': [0, 1, 10, 11],     # List of vertex IDs
    'edges': [(0,1), (1,11), ...], # List of edge tuples  
    'weight': 4,                    # Number of edges
    'id': 0                        # Unique identifier
}
```

### Cluster Representation  
```python
cluster = {
    'loop_ids': [0, 1, 2],              # IDs of loops in cluster
    'multiplicities': {0: 2, 1: 1, 2: 1}, # How many times each loop appears
    'weight': 7,                         # Total weight (sum of loop weights × multiplicities)
    'total_loops': 4                     # Total number of loop instances
}
```

### Interaction Graph
```python
interaction_graph = {
    0: [1, 3, 5],    # Loop 0 is incompatible with loops 1, 3, 5
    1: [0, 2, 4],    # Loop 1 is incompatible with loops 0, 2, 4
    ...
}
```

## 🔄 Typical Usage Patterns

### 1. Find All Loops on a Site
```python
from optimized_loop_enumeration import OptimizedLoopEnumerator, create_periodic_square_lattice

L = 10
site = 0  
max_weight = 7

adj_matrix = create_periodic_square_lattice(L)
enumerator = OptimizedLoopEnumerator(adj_matrix)
loops = enumerator.find_loops_supported_on_vertex(site, max_weight)

print(f"Found {len(loops)} loops")
```

### 2. Find All Connected Clusters
```python
from cluster_enumeration_complete import ClusterEnumerator, create_periodic_square_lattice

L = 10
site = 0
max_weight = 7

adj_matrix = create_periodic_square_lattice(L)
enumerator = ClusterEnumerator(adj_matrix)
clusters = enumerator.enumerate_connected_clusters(site, max_weight)

print(f"Found {len(clusters)} clusters")
```

### 3. Analyze Results
```python
# Group loops by weight
from collections import defaultdict

by_weight = defaultdict(list)
for loop in loops:
    by_weight[loop['weight']].append(loop)

for weight, loops_of_weight in by_weight.items():
    print(f"Weight {weight}: {len(loops_of_weight)} loops")
```

## ⚠️ Important Notes

1. **Coordinate System:** Vertices are numbered 0 to L²-1 in row-major order
2. **Edges:** Always stored as tuples (u,v) with u < v  
3. **Weight:** Always refers to number of edges, not vertices
4. **Periodic Boundaries:** All lattices use periodic boundary conditions (torus topology)
5. **Degree Constraint:** Loops require ALL vertices to have degree ≥ 2 in the subgraph