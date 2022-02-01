import graph_tool.stats as gts

def edge_weights(g, hub_penalty, inverse=False):
    avdeg = gts.vertex_average(g, "total")[0]
    weights = g.new_edge_property("double", val=avdeg)
    if hub_penalty <= 0:
        return weights
    if hub_penalty > 1:
        raise ValueError("Invalid hub penalty {}.".format(hub_penalty))
    for e in g.edges():
        edge_avdeg = float(e.source().out_degree() + e.target().out_degree()) / 2.0
        penalized_weight = (1.0 - hub_penalty) * avdeg + hub_penalty * edge_avdeg
        if inverse:
            weights[e] = 1.0 / penalized_weight
        else:
            weights[e] = penalized_weight
    return weights