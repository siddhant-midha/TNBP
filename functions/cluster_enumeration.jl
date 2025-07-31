#=
This file contains Python code for enumerating loops supported on a given vertex.
The actual Python implementation is in loop_enumeration.py and test_loop_enumeration.py

Usage:
    python test_loop_enumeration.py

Description:
A loop is defined as a connected subgraph where each vertex has degree >= 2.
The implementation uses depth-first search to find all such loops starting from
a given vertex with weight (number of edges) <= m.

The LoopEnumerator class provides:
- find_loops_supported_on_vertex(vertex, max_weight): Find all loops containing the vertex
- Support for periodic square lattices
- Loop deduplication and validation
- DFS-based enumeration with edge tracking

Test coverage includes:
- 3x3 and 4x4 periodic square lattices
- Different vertex types (corner, edge, center)
- Various weight limits
- Performance testing on larger lattices
- Loop connectivity verification

To run tests:
    cd /path/to/TNBP
    python test_loop_enumeration.py
=#