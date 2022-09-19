def make_node_id_map(g):
    mapping = {}
    for node in range(g.num_vertices()):
        mapping[g.vertex_properties['internal_id'][node]] = node
    return mapping

def add_edges(g, edge_list):
    """
    edge list is [{"fom":..., "to":...}, ...]
    """
    mapping = make_node_id_map(g)
    edge_id_list = []
    for edge in edge_list:
        a = mapping[edge['from']]
        b = mapping[edge['to']]
        edge_id_list.append((a, b))
    g.add_edge_list(edge_list)
    return g