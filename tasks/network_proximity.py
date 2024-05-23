from tasks.task_hook import TaskHook
from tasks.util.custom_network import add_edges, remove_ppi_edges, filter_proteins
from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.edge_weights import edge_weights
import os.path
import graph_tool as gt
import graph_tool.topology as gtt
import sys
import numpy as np


def network_proximity(task_hook: TaskHook):

    # Type: List of str
    # Semantics: Names of the seed proteins. Use UNIPROT AC for host proteins, and
    #            names of the format SARS_CoV2_<IDENTIFIER> (tree_edge.g., SARS_CoV2_ORF6) for
    #            virus proteins.
    # Reasonable default: None, has to be selected by user via "select for analysis"
    #            utility in frontend.
    # Acceptable values: UNIPROT ACs, identifiers of viral proteins.
    seeds = task_hook.parameters["seeds"]
    nodes_not_in_lcc = {'D3W0D1', 'O15178', 'O60542', 'O60609', 'O60882', 'Q5T4W7', 'Q6UVW9', 'Q6UW32', 'Q6UWQ7',
                        'Q6UXB1', 'Q7RTX0', 'Q7RTX1', 'Q8TE23', 'Q9GZZ7', 'Q9H2W2', 'Q9H665', 'Q9NZW4'}
    # seeds = [node for node in set(seeds).difference(nodes_not_in_lcc)]

    # Type: bool
    # Semantics: Sepcifies whether should be included in the analysis when ranking drugs.
    # Example: False.
    # Reasonable default: False.
    # Has no effect unless trust_rank.py is used for ranking drugs.
    include_non_approved_drugs = task_hook.parameters.get("include_non_approved_drugs", False)

    # Type: int.
    # Semantics: Number of random seed sets for computing Z-scores.
    # Example: 32.
    # Reasonable default: 32.
    # Acceptable values: Positive integers.
    num_random_seed_sets = task_hook.parameters.get("num_random_seed_sets", 32)

    # Type: int.
    # Semantics: Number of random drug target sets for computing Z-scores.
    # Example: 32.
    # Reasonable default: 32.
    # Acceptable values: Positive integers.
    num_random_drug_target_sets = task_hook.parameters.get("num_random_drug_target_sets", 32)

    # Type: int.
    # Semantics: Number of returned drugs.
    # Example: 20.
    # Reasonable default: 20.
    # Acceptable values: integers n with n > 0.
    result_size = task_hook.parameters.get("result_size", 20)

    # Type: int.
    # Semantics: All nodes with degree > max_deg * g.num_vertices() are ignored.
    # Example: 39.
    # Reasonable default: sys.maxsize.
    # Acceptable values: Positive integers.
    max_deg = task_hook.parameters.get("max_deg", sys.maxsize)

    # Type: float.
    # Semantics: Penalty parameter for hubs. Set edge weight to 1 + (hub_penalty / 2) (e.source.degree + e.target.degree)
    # Example: 0.5.
    # Reasonable default: 0.
    # Acceptable values: Floats between 0 and 1.
    hub_penalty = task_hook.parameters.get("hub_penalty", 0.0)

    # Type: int.
    # Semantics: Number of threads used for running the analysis.
    # Example: 1.
    # Reasonable default: 1.
    # Note: We probably do not want to expose this parameter to the user.
    num_threads = task_hook.parameters.get("num_threads", 1)
    
    ppi_dataset = task_hook.parameters.get("ppi_dataset")

    pdi_dataset = task_hook.parameters.get("pdi_dataset")

    search_target = task_hook.parameters.get("target", "drug-target")

    filter_paths = task_hook.parameters.get("filter_paths", True)

    custom_edges = task_hook.parameters.get("custom_edges", False)
    
    no_default_edges =    no_default_edges = task_hook.parameters.get("exclude_drugstone_ppi_edges", False)
    
    custom_nodes = task_hook.parameters.get("network_nodes", False)

    node_name_attribute = "internal_id"  # nodes in the input network which is created from RepoTrialDB have primaryDomainId as name attribute
    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)

    # Parsing input file.
    task_hook.set_progress(0.0 / 8, "Parsing input.")

    id_space = task_hook.parameters["config"].get("identifier", "symbol")

    filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset['licenced'] or pdi_dataset['licenced']:
        filename += "_licenced"
    filename = os.path.join(task_hook.data_directory, filename + ".gt")
    # g, seed_ids, _, drug_ids = read_graph_tool_graph(file_path, seeds, "", "", max_deg, False, True, include_non_approved_drugs)
    g, seed_ids, drug_ids = read_graph_tool_graph(filename, seeds, id_space, max_deg, True, include_non_approved_drugs, target=search_target)
    
    if custom_edges:
      if no_default_edges:
        # clear all edges with type "protein-protein"
        g = remove_ppi_edges(g)
      g = add_edges(g, custom_edges)
    
    if custom_nodes:
      # remove all nodes with internal_id not in custom_nodes from g
      g, seed_ids, drug_ids = filter_proteins(g, custom_nodes, drug_ids, seeds)
    
    # Computing edge weights.
    task_hook.set_progress(1.0 / 8, "Computing edge weights.")
    weights = edge_weights(g, hub_penalty)

    # Delete drug targets not in LCC.
    task_hook.set_progress(2.0 / 8, "Deleting drug targets not in LCC.")
    drug_targets = {drug_id: g.get_all_neighbors(drug_id) for drug_id in drug_ids}
    for drug_id in drug_ids:
        deleted_targets = []
        targets = drug_targets[drug_id]
        for i in range(len(targets)):
            # if g.vertex_properties["name"][targets[i]] in nodes_not_in_lcc:
            if g.vertex_properties[node_name_attribute][targets[i]] in nodes_not_in_lcc:
                deleted_targets.append(i)
        drug_targets[drug_id] = np.delete(targets, deleted_targets)

    # Compute all shortest path distances.
    task_hook.set_progress(3.0 / 8, "Computing all shortest path distances.")
    distances = None
    if hub_penalty == 0:
        # distances = np.load(os.path.join(task_hook.data_directory, "host-host-distances.npy"))  # need to create this file and also rename it to protein-protein-distances.npy
        distances = gtt.shortest_distance(g)
    else:
        distances = gtt.shortest_distance(g, weights=weights)

    # Compute network proximities.
    task_hook.set_progress(4.0 / 8, "Computing network proximities.")
    proximities = {drug_id : 0 for drug_id in drug_ids}
    for drug_id in drug_ids:
        distance = 0.0
        for drug_target in drug_targets[drug_id]:
            distance += min([distances[drug_target][seed_id] for seed_id in seed_ids])
        if len(drug_targets[drug_id]) == 0:
            proximities[drug_id] = np.inf
        else:
            proximities[drug_id] = distance / float(len(drug_targets[drug_id]))

    # Compute background distribution.
    task_hook.set_progress(5.0 / 8, "Computing background distribution")
    min_num_targets = min([len(drug_targets[drug_id]) for drug_id in drug_ids])
    max_num_targets = max([len(drug_targets[drug_id]) for drug_id in drug_ids])
    node_ids_in_lcc = [node for node in range(g.num_vertices()) if
                       # not g.vertex_properties["name"][node] in nodes_not_in_lcc]
                       not g.vertex_properties[node_name_attribute][node] in nodes_not_in_lcc]
    background_distribution = []
    num_seeds = len(seed_ids)
    for i in range(num_random_seed_sets):
        np.random.shuffle(node_ids_in_lcc)
        random_seed_ids = node_ids_in_lcc[:num_seeds]
        for k in range(num_random_drug_target_sets):
            np.random.shuffle(node_ids_in_lcc)
            random_drug_targets = node_ids_in_lcc[:np.random.randint(min_num_targets, max_num_targets + 1)]
            distance = 0.0
            for drug_target in random_drug_targets:
                distance += min([distances[drug_target][seed_id] for seed_id in random_seed_ids])
            background_distribution.append(distance / float(len(random_drug_targets)))
    background_mean = np.mean(background_distribution)
    background_std = np.std(background_distribution)

    # Apply Z-score transformation.
    task_hook.set_progress(6.0 / 8, "Applying Z-score transformation.")
    drugs_with_z_scores = [(drug_id, (proximities[drug_id] - background_mean) / background_std) for drug_id in drug_ids]

    task_hook.set_progress(7.0 / 8, "Formatting results.")
    best_drugs = [item for item in sorted(drugs_with_z_scores, key=lambda item: item[1])[:result_size]]
    best_drugs_ids = [item[0] for item in best_drugs]
    seed_ids = list(set(seed_ids))
        # Concatenate best result candidates with seeds and compute induced subgraph.
    # since the result size filters out nodes, the result network is not complete anymore.
    # Therefore, it is necessary to find the shortest paths to the found nodes in case intermediate nodes have been removed. 
    # This leads to more nodes in the result that given in the threshold, but the added intermediate nodes help to understand the output 
    # and are marked as intermediate nodes.
    intermediate_nodes = set()

    returned_edges = set()
    returned_nodes = set(seed_ids) # return seed_ids in any case

    # return only the path to a drug with the shortest distance
    if filter_paths:
        for candidate in best_drugs_ids:
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
        for candidate in best_drugs_ids:
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
    subgraph = {
        "nodes": [g.vertex_properties[node_name_attribute][node] for node in returned_nodes],
        "edges": [{"from": g.vertex_properties[node_name_attribute][source], "to": g.vertex_properties[node_name_attribute][target]} for source, target in returned_edges],
        }

    # Compute node attributes.
    node_types = {g.vertex_properties[node_name_attribute][node]: g.vertex_properties["type"][node] for node in returned_nodes}
    is_seed = {g.vertex_properties[node_name_attribute][node]: node in set(seed_ids) for node in returned_nodes}
    returned_scores = {g.vertex_properties[node_name_attribute][node]: None for node in returned_nodes}

    for node, score in best_drugs:
        # returned_scores[g.vertex_properties["name"][node]] = score
        returned_scores[g.vertex_properties[node_name_attribute][node]] = score

    # accepted_candidates are needed to comply with the output format of "scores_to_results"
    accepted_candidates = [x for x in subgraph['nodes'] if x[:2] == 'dr']
    
    task_hook.set_results({
        "network": subgraph,
        'intermediate_nodes': list(intermediate_nodes),
        'target_nodes': accepted_candidates,
        "node_attributes":
            {
                "node_types": node_types,
                "is_seed": is_seed,
                "scores": returned_scores
            },
        'gene_interaction_dataset': ppi_dataset,
        'drug_interaction_dataset': pdi_dataset,
    })
