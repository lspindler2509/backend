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
        a = mapping[edge['from']] if edge['from'] in mapping else False
        b = mapping[edge['to']] if edge['to'] in mapping else False
        if a and b:
            edge_id_list.append((a, b, 'protein-protein'))
    e_type = g.edge_properties["type"]
    g.add_edge_list(edge_id_list, eprops=[e_type])
    return g