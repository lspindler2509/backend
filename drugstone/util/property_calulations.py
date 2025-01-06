import networkx as nx

def name2index(g, node_name_attribute="internal_id"):
    """
    Create a mapping from gene name to vertex index.
    """
    index2name = g.vertex_properties[node_name_attribute]
    return {index2name[v]: v for v in g.iter_vertices()}

def find_vertices(ids, g ,mapping):
    """Find vertices in the graph for given IDs."""
    found_vertices = {}
    for node in ids:
        if node.startswith('dr'):
            continue
        vertex = mapping.get(node, None)
        if vertex:
            found_vertices[node] = g.vertex(vertex)
        else:
            found_vertices[node] = None

    return found_vertices

def build_nx_graph(edges, valid_ids):
    """Build a NetworkX graph from valid IDs and edges."""
    nx_graph = nx.Graph()
    nx_graph.add_nodes_from(valid_ids)
    for edge in edges:
        source = edge['from']
        target = edge['to']
        if source in valid_ids and target in valid_ids:
            nx_graph.add_edge(source, target)
    return nx_graph

def calculate_network_properties(nx_graph, node_id, degree_in_ppi):
    """Calculate network properties for a single node."""
    nx_degree = nx_graph.degree[node_id]
    nx_clustering = nx.clustering(nx_graph, node_id)
    spd = nx_degree / degree_in_ppi if degree_in_ppi > 0 else 0
    return nx_degree, nx_clustering, spd

def calculate_properties_id_based(ids, g, edges, calculateProperties = True):
    if not calculateProperties:
        return {node: {} for node in ids}

    if not g:
        print("No graph given")
        return {}

    properties = {}
    mapping = name2index(g)
    found_vertices = find_vertices(ids, g ,mapping)

    valid_ids = {id for id, vertex in found_vertices.items() if vertex}
    nx_graph = build_nx_graph(edges, valid_ids)

    for node in ids:
        if node.startswith('dr'):
            continue
        properties[node] = {}
        if node in valid_ids:
            vertex = found_vertices.get(node)
            degree_in_ppi = calculate_filtered_degree(g, vertex, "protein-protein") if vertex else 0
            properties[node]['degree_in_ppi'] = degree_in_ppi
            nx_degree, nx_clustering, spd = calculate_network_properties(nx_graph, node, degree_in_ppi)
            properties[node]['degree_in_network'] = nx_degree
            properties[node]['local_clustering_coefficient'] = nx_clustering
            properties[node]['SPD'] = spd
        else:
            print(f"Skipping node ID {node} as it is not in the graph.")
    
    return properties

def calculate_properties(nodes, g, identifier, edges, calculateProperties = True):
    if not calculateProperties:
        for node in nodes:
            node.setdefault('properties', {})
        return nodes

    if not g:
        print("No graph given")
        return nodes

    ids = [node[identifier][0] for node in nodes if identifier in node]
    mapping = name2index(g)
    found_vertices = find_vertices(ids, g ,mapping)

    valid_ids = {id for id, vertex in found_vertices.items() if vertex}
    nx_graph = build_nx_graph(edges, valid_ids)
    for node in nodes:
        node.setdefault('properties', {})
        id = node[identifier][0] if identifier in node else None
        if id and id in valid_ids:
            vertex = found_vertices.get(id)
            degree_in_ppi = calculate_filtered_degree(g, vertex, "protein-protein")
            node['properties']['degree_in_ppi'] = degree_in_ppi
            nx_degree, nx_clustering, spd = calculate_network_properties(nx_graph, id, degree_in_ppi)
            node['properties']['degree_in_network'] = nx_degree
            node['properties']['local_clustering_coefficient'] = nx_clustering
            node['properties']['SPD'] = spd
        else:
            print(f"Skipping node ID {id} as it is not in the graph.")
    return nodes



def calculate_filtered_degree(g, vertex, target_type):
    """
    Calculates the degree of a vertex considering only edges of a specific type.
    """
    edge_type = g.edge_properties["type"]
    # Count edges connected to the vertex that match the target type
    return sum(1 for edge in vertex.all_edges() if edge_type[edge] == target_type)
