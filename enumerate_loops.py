import networkx as nx
import itertools
from collections import deque

def canonical_edge(u, v):
    """
    Creates a canonical representation of an edge by sorting the nodes.
    This helps in comparing edge sets regardless of node order.
    
    Args:
        u: The first node of the edge.
        v: The second node of the edge.
        
    Returns:
        A tuple (node1, node2) where node1 <= node2.
    """
    return tuple(sorted((u, v)))

def find_simple_cycles_through_node(G, start_node, max_weight):
    """
    Finds all simple cycles containing start_node with weight (edge count)
    less than or equal to max_weight.

    This function uses a recursive backtracking algorithm (a modified DFS).

    Args:
        G (nx.Graph): The input graph.
        start_node: The node that all cycles must pass through.
        max_weight (int): The maximum number of edges the cycles can have.

    Returns:
        A list of frozensets, where each frozenset represents the
        canonical edges of a found cycle.
    """
    cycles = set()
    
    def find_cycles_recursive(current_node, path, path_set):
        if len(path) > max_weight:
            return

        for neighbor in G.neighbors(current_node):
            if neighbor == start_node:
                if len(path) > 2:
                    edge_set = set()
                    for i in range(len(path) - 1):
                        edge_set.add(canonical_edge(path[i], path[i+1]))
                    edge_set.add(canonical_edge(path[-1], start_node))
                    
                    if len(edge_set) <= max_weight:
                        cycles.add(frozenset(edge_set))
            elif neighbor not in path_set:
                path.append(neighbor)
                path_set.add(neighbor)
                find_cycles_recursive(neighbor, path, path_set)
                path_set.remove(neighbor)
                path.pop()

    find_cycles_recursive(start_node, [start_node], {start_node})
    return list(cycles)

def find_all_loops(G, start_node, max_weight):
    """
    Finds all "loops" supported on a given vertex with a maximum weight.

    A loop is a connected subgraph where every vertex has a degree of at least 2.

    Args:
        G (nx.Graph): The input graph from networkx.
        start_node: The vertex that all loops must contain.
        max_weight (int): The maximum allowed weight for a loop.

    Returns:
        A list of networkx.Graph objects, each representing a valid loop.
    """
    if start_node not in G:
        raise ValueError("The specified start_node is not in the graph.")

    found_loops_edges = set()
    process_queue = deque()

    initial_cycles = find_simple_cycles_through_node(G, start_node, max_weight)
    for cycle_edges in initial_cycles:
        if cycle_edges not in found_loops_edges:
            found_loops_edges.add(cycle_edges)
            process_queue.append(cycle_edges)
            
    while process_queue:
        current_loop_edges = process_queue.popleft()
        current_loop_graph = nx.Graph(list(current_loop_edges))
        current_loop_nodes = set(current_loop_graph.nodes())
        
        for u, v in itertools.combinations(current_loop_nodes, 2):
            nodes_to_avoid = current_loop_nodes - {u, v}
            G_temp = G.copy()
            G_temp.remove_nodes_from(nodes_to_avoid)
            
            max_path_len = max_weight - len(current_loop_edges)
            if max_path_len < 1:
                continue

            paths = nx.all_simple_paths(G_temp, source=u, target=v, cutoff=max_path_len + 1)
            for path in paths:
                if len(path) > 1:
                    new_path_edges = {canonical_edge(path[i], path[i+1]) for i in range(len(path) - 1)}
                    new_loop_edges = current_loop_edges.union(new_path_edges)

                    if len(new_loop_edges) <= max_weight and new_loop_edges not in found_loops_edges:
                        found_loops_edges.add(new_loop_edges)
                        process_queue.append(new_loop_edges)

    result_graphs = [nx.Graph(list(edge_set)) for edge_set in found_loops_edges]
    return result_graphs

def _generate_all_loops_in_graph(G, max_weight):
    """
    Helper function to generate a database of all possible loops in the entire
    graph up to a given max_weight. This is a computationally intensive step.

    Returns:
        A dictionary (loop_db) mapping a loop's edge set (frozenset) to its
        properties (weight, nodes, graph object).
    """
    # Step 1: Find all simple cycles by iterating over every node as a start node.
    # Using a set handles duplicates automatically.
    print("Step 1/3: Finding all simple cycles in the graph (this may take a while)...")
    all_simple_cycles = set()
    for i, node in enumerate(G.nodes()):
        # Provide some progress feedback for large graphs
        if i % 10 == 0:
            print(f"  - Scanning for cycles from node {i+1}/{len(G.nodes())}...")
        cycles_through_node = find_simple_cycles_through_node(G, node, max_weight)
        all_simple_cycles.update(cycles_through_node)
    
    # Step 2: Expand simple cycles into all complex loops using the "ear-addition" method.
    print(f"Step 2/3: Found {len(all_simple_cycles)} simple cycles. Expanding to complex loops...")
    found_loops_edges = set(all_simple_cycles)
    process_queue = deque(all_simple_cycles)
    
    # This is the same expansion logic from find_all_loops
    while process_queue:
        current_loop_edges = process_queue.popleft()
        current_loop_graph = nx.Graph(list(current_loop_edges))
        current_loop_nodes = set(current_loop_graph.nodes())
        
        for u, v in itertools.combinations(current_loop_nodes, 2):
            nodes_to_avoid = current_loop_nodes - {u, v}
            G_temp = G.copy()
            G_temp.remove_nodes_from(nodes_to_avoid)
            
            max_path_len = max_weight - len(current_loop_edges)
            if max_path_len < 1:
                continue
            
            paths = nx.all_simple_paths(G_temp, source=u, target=v, cutoff=max_path_len + 1)
            for path in paths:
                if len(path) > 1:
                    new_path_edges = {canonical_edge(path[i], path[i+1]) for i in range(len(path) - 1)}
                    new_loop_edges = current_loop_edges.union(new_path_edges)

                    if len(new_loop_edges) <= max_weight and new_loop_edges not in found_loops_edges:
                        found_loops_edges.add(new_loop_edges)
                        process_queue.append(new_loop_edges)

    # Step 3: Build a database for quick lookup of loop properties.
    print(f"Step 3/3: Generated {len(found_loops_edges)} unique loops. Building database...")
    loop_db = {}
    for edge_set in found_loops_edges:
        loop_graph = nx.Graph(list(edge_set))
        loop_db[edge_set] = {
            'weight': loop_graph.number_of_edges(),
            'nodes': set(loop_graph.nodes()),
            'graph': loop_graph
        }
    print("Loop database created.\n")
    return loop_db

def find_all_connected_clusters(G, start_node, max_weight):
    """
    Enumerates all connected clusters supported on a given site with a total
    weight at most max_weight.

    Args:
        G (nx.Graph): The input graph.
        start_node: The site the clusters must be supported on.
        max_weight (int): The maximum total weight of a cluster.

    Returns:
        A list of clusters. Each cluster is a list of networkx.Graph objects.
    """
    # Step 1: Generate the universe of all possible loops up to max_weight
    loop_db = _generate_all_loops_in_graph(G, max_weight)

    # Step 2: Identify loops supported on the start_node to seed the search
    initial_loop_ids = {
        lid for lid, data in loop_db.items() 
        if start_node in data['nodes']
    }

    # Step 3: Use a BFS-style search to find all valid connected clusters
    print("Searching for connected clusters...")
    found_clusters = set()  # Stores canonical representations of clusters to avoid duplicates
    queue = deque()

    # Define a deterministic key for sorting frozensets of edges
    frozenset_sort_key = lambda fs: tuple(sorted(list(fs)))

    # Initialize queue with single-loop clusters from the initial set
    for lid in initial_loop_ids:
        # A cluster is represented by a tuple of loop IDs (frozensets).
        # We sort it to create a canonical representation for the set.
        cluster_tuple = (lid,)
        weight = loop_db[lid]['weight']
        if cluster_tuple not in found_clusters:
            found_clusters.add(cluster_tuple)
            queue.append((cluster_tuple, weight))

    # BFS search to grow clusters
    processed_count = 0
    while queue:
        current_cluster, current_weight = queue.popleft()
        processed_count += 1
        if processed_count % 100 == 0:
            print(f"  - Processed {processed_count} potential clusters...")
            
        nodes_in_cluster = set().union(*(loop_db[lid]['nodes'] for lid in current_cluster))

        # Try to add every possible loop from our database to the current cluster
        for cand_id, cand_data in loop_db.items():
            new_weight = current_weight + cand_data['weight']
            if new_weight > max_weight:
                continue

            # Check for incompatibility (must share at least one node)
            # This ensures the interaction graph of the new cluster remains connected.
            if not nodes_in_cluster.intersection(cand_data['nodes']):
                continue

            # Form the new cluster and its canonical representation
            new_cluster_list = list(current_cluster) + [cand_id]
            new_cluster_canonical = tuple(sorted(new_cluster_list, key=frozenset_sort_key))

            if new_cluster_canonical in found_clusters:
                continue

            found_clusters.add(new_cluster_canonical)
            queue.append((new_cluster_canonical, new_weight))

    # Final Step: Convert results from loop IDs to actual graph objects
    final_clusters = []
    for cluster_tuple in found_clusters:
        cluster_graphs = [loop_db[lid]['graph'] for lid in cluster_tuple]
        final_clusters.append(cluster_graphs)

    return final_clusters

# find all connected clusters supported on any vertex in the graph
# and with a total weight at most max_weight.
def find_connected_clusters(G, max_weight):

    all_clusters = []
    for node in G.nodes():
        clusters = find_all_connected_clusters(G, node, max_weight)
        # check redundancy
        all_clusters.append(clusters)

    return all_clusters

def pretty_print_loop(loop_graph, indent="  "):
    """Helper function to print the properties of a loop in a readable format."""
    weight = loop_graph.number_of_edges()
    nodes = sorted(list(loop_graph.nodes()))
    edges = sorted([tuple(sorted(e)) for e in loop_graph.edges()])
    print(f"{indent}- Loop with weight {weight}:")
    print(f"{indent}  Nodes ({len(nodes)}): {nodes}")
    print(f"{indent}  Edges ({len(edges)}): {edges}")
    
def pretty_print_cluster(cluster):
    """Helper function to print a cluster and its properties."""
    total_weight = sum(loop.number_of_edges() for loop in cluster)
    print(f"- Cluster of size {len(cluster)} | Total Weight: {total_weight}")
    for i, loop_graph in enumerate(cluster):
        print(f"  Loop #{i+1}:")
        pretty_print_loop(loop_graph, indent="    ")

# --- Main block to test the code ---
if __name__ == '__main__':
    # Using a smaller lattice for cluster enumeration as it's more intensive
    width, height = 10, 10
    G = nx.grid_2d_graph(width, height, periodic=True)
    print(f"Created a {width}x{height} periodic square lattice.")
    print("-" * 40)

    # Set the parameters
    start_vertex = (0, 0)
    # A modest max weight is recommended as the number of clusters grows very rapidly.
    max_cluster_weight = 8
    
    print(f"Finding all connected clusters supported on vertex {start_vertex} with total weight <= {max_cluster_weight}\n")

    try:
        # Run the main cluster-finding function
        found_clusters = find_all_connected_clusters(G, start_vertex, max_cluster_weight)

        # Print the results
        print("-" * 40)
        print(f"\nFound {len(found_clusters)} unique connected clusters.")
        
        # Sort clusters by total weight, then by size for cleaner output
        found_clusters.sort(key=lambda c: (sum(l.number_of_edges() for l in c), len(c)))
        
        for i, cluster in enumerate(found_clusters):
            print("\n" + "=" * 20 + f" Cluster #{i+1} " + "=" * 20)
            pretty_print_cluster(cluster)
            
    except ValueError as e:
        print(f"Error: {e}")
