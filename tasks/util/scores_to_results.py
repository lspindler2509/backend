import graph_tool.topology as gtt


def scores_to_results(
        target,
        result_size,
        g,
        seed_ids,
        drug_ids,
        scores,
        ppi_dataset,
        pdi_dataset,
        filterPaths
):

    r"""Transforms the scores to the required result format."""

    node_name_attribute = "drugstone_id"  # nodes in the input network which is created from RepoTrialDB have primaryDomainId as name attribute
    candidates = []
    # if strain_or_drugs == "drugs":
    if target == "drug":
        candidates = [(node, scores[node]) for node in drug_ids if scores[node] > 0]
    else:
        candidates = [(node, scores[node]) for node in range(g.num_vertices()) if scores[node] > 0 and node not in set(seed_ids)]
    best_candidates = [item[0] for item in sorted(candidates, key=lambda item: item[1], reverse=True)[:result_size]]
    print(f'Candidate list length: {len(best_candidates)}')

    # Concatenate best result candidates with seeds and compute induced subgraph.
    # since the result size filters out nodes, the result network is not complete anymore.
    # Therefore, it is necessary to find the shortest paths to the found nodes in case intermediate nodes have been removed. 
    # This leads to more nodes in the result that given in the threshold, but the added intermediate nodes help to understand the output 
    # and are marked as intermediate nodes.
    intermediate_nodes = set()

    returned_edges = set()
    returned_nodes = set(seed_ids) # return seed_ids in any case

    # return only the path to a drug with the shortest distance
    if filterPaths:
        for candidate in best_candidates:
            distances = gtt.shortest_distance(g, candidate, seed_ids)
            closest_distance_mean = sum(distances) / len(distances)

            for index, seed_id in enumerate(seed_ids):
                if distances[index] > closest_distance_mean:
                    continue
                vertices, edges = gtt.shortest_path(g, candidate, seed_id)

                drug_in_path = False
                for vertex in vertices:
                    if g.vertex_properties["type"][int(vertex)] == "Drug" and vertex != candidate:
                        drug_in_path = True
                        break
                if drug_in_path:
                    continue

                for vertex in vertices:
                    if int(vertex) not in returned_nodes:
                        # inserting intermediate node in order to make result comprehensive
                        intermediate_nodes.add(g.vertex_properties[node_name_attribute][int(vertex)])
                        returned_nodes.add(int(vertex))
                for edge in edges:
                    if ((edge.source(), edge.target()) not in returned_edges) or ((edge.target(), edge.source()) not in returned_edges):
                        returned_edges.add((edge.source(), edge.target()))
    else:
        for candidate in best_candidates:
            for index, seed_id in enumerate(seed_ids):
                vertices, edges = gtt.shortest_path(g, candidate, seed_id)

                drug_in_path = False
                for vertex in vertices:
                    if g.vertex_properties["type"][int(vertex)] == "Drug" and vertex != candidate:
                        drug_in_path = True
                        break
                if drug_in_path:
                    continue

                for vertex in vertices:
                    if int(vertex) not in returned_nodes:
                        # inserting intermediate node in order to make result comprehensive
                        intermediate_nodes.add(g.vertex_properties[node_name_attribute][int(vertex)])
                        returned_nodes.add(int(vertex))
                for edge in edges:
                    if ((edge.source(), edge.target()) not in returned_edges) or ((edge.target(), edge.source()) not in returned_edges):
                        returned_edges.add((edge.source(), edge.target()))
    print(f'Returned nodes number: {len(returned_nodes)}')
    subgraph = {
        "nodes": [g.vertex_properties[node_name_attribute][node] for node in returned_nodes],
        "edges": [{"from": g.vertex_properties[node_name_attribute][source], "to": g.vertex_properties[node_name_attribute][target]} for source, target in returned_edges],
        }

    # Compute node attributes.
    node_types = {g.vertex_properties[node_name_attribute][node]: g.vertex_properties["type"][node] for node in returned_nodes}
    is_seed = {g.vertex_properties[node_name_attribute][node]: node in set(seed_ids) for node in returned_nodes}
    returned_scores = {g.vertex_properties[node_name_attribute][node]: scores[node] for node in returned_nodes}

    return {
        "network": subgraph,
        'intermediate_nodes': list(intermediate_nodes),
        "node_attributes":
            {
                "node_types": node_types,
                "is_seed": is_seed,
                "scores": returned_scores
            },
        'gene_interaction_dataset': ppi_dataset,
        'drug_interaction_dataset': pdi_dataset,
    }
