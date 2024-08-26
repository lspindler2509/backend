from tasks.task_hook import TaskHook
from tasks.util.custom_network import add_edges, remove_ppi_edges, filter_proteins
from tasks.util.steiner_tree import steiner_tree
from tasks.util.find_bridges import find_bridges
from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.edge_weights import edge_weights
import os.path
import graph_tool as gt
import graph_tool.util as gtu
import sys


def multi_steiner(task_hook: TaskHook):
    
    # Type: List of str
    # Semantics: Names of the seed proteins. Use UNIPROT AC for host proteins, and
    #            names of the format SARS_CoV2_<IDENTIFIER> (tree_edge.g., SARS_CoV2_ORF6) for
    #            virus proteins.
    # Reasonable default: None, has to be selected by user via "select for analysis"
    #            utility in frontend.
    # Acceptable values: UNIPROT ACs, identifiers of viral proteins.
    seeds = task_hook.parameters["seeds"]
    # Since we get the LCC of the graphs from read_graph_tool_graph, no need for the following two lines
    # nodes_not_in_lcc = {'D3W0D1', 'O15178', 'O60542', 'O60609', 'O60882', 'Q5T4W7', 'Q6UVW9', 'Q6UW32', 'Q6UWQ7', 'Q6UXB1', 'Q7RTX0', 'Q7RTX1', 'Q8TE23', 'Q9GZZ7', 'Q9H2W2', 'Q9H665', 'Q9NZW4'}
    # seeds = [node for node in set(seeds).difference(nodes_not_in_lcc)]
    seeds.sort()

    # Type: str.
    # Semantics: The virus strain for which the analysis should be run.
    # Example: "SARS_CoV2"
    # Reasonable default: None, has to be specified by the caller.
    ### Acceptable values: "SARS_CoV2", ...
    # strain_or_drugs = task_hook.parameters.get("strain_or_drugs", "SARS_CoV2")
    # Acceptable values: "PPI", "PPDr"
    # target_or_drugs = task_hook.parameters.get("target_or_drugs", "PPI")

    # Type: list of str.
    # Semantics: The datasets which should be considered for the analysis.
    # Example: ["Krogan", "TUM"].
    # Note: If empty, all available datasets are used.
    # Reasonable default: [].
    # Acceptable values: "Krogan", "TUM".
    # datasets = task_hook.parameters.get("datasets", [])

    # Type: list of str.
    # Semantics: Virus-host g_edge types which should be ignored for the analysis.
    # Example: ["Overexpression"].
    # Note: If empty, all available g_edge types are used.
    # Reasonable default: [].
    # Acceptable values: "AP-MS", "overexpression".
    # ignored_edge_types = task_hook.parameters.get("ignored_edge_types", [])
    
    # Type bool.
    # Semantics: Ignore viral proteins which are not selected as seeds.
    # Example: False.
    # Reasonable default: False.
    # Has no effect when the algorithm is used for ranking drugs.
    # ignore_non_seed_baits = task_hook.parameters.get("ignore_non_seed_baits", False)

    # Type: int.
    # Semantics: Number of steiner trees to return.
    # Example: 10.
    # Reasonable default: 5.
    # Acceptable values: integers n with 0 < n < 25.
    num_trees = task_hook.parameters.get("num_trees", 5)
    
    # Type: float.
    # Semantics: The error tolerance of the subsequent Steiner trees w.r.t. the first one in percent.
    # Example: 10
    # Reasonable default: 10.
    # Acceptable values: floats x >= 0. 
    tolerance = task_hook.parameters.get("tolerance", 10)
    
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

    node_name_attribute = "internal_id" # nodes in the input network which is created from RepoTrialDB have primaryDomainId as name attribute

    custom_edges = task_hook.parameters.get("custom_edges", False)
    
    no_default_edges =    no_default_edges = task_hook.parameters.get("exclude_drugstone_ppi_edges", False)
    
    custom_nodes = task_hook.parameters.get("network_nodes", False)

    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    # Parsing input file.
    task_hook.set_progress(0 / (float(num_trees + 3)), "Parsing input.")

    id_space = task_hook.parameters["config"].get("identifier", "symbol")

    filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset['licenced'] or pdi_dataset['licenced']:
        filename += "_licenced"
    filename = os.path.join(task_hook.data_directory, filename + ".gt")
    g, seed_ids, _ = read_graph_tool_graph(filename, seeds, id_space, max_deg, target=search_target)

    if custom_edges:
      if no_default_edges:
        # clear all edges with type "protein-protein"
        g = remove_ppi_edges(g)
      g = add_edges(g, custom_edges)
    
    if custom_nodes:
      # remove all nodes with internal_id not in custom_nodes from g
      g, seed_ids, drug_ids = filter_proteins(g, custom_nodes, drug_ids, seeds)

    seed_map = {g.vertex_properties[node_name_attribute][node]: node for node in seed_ids}
    task_hook.set_progress(1 / (float(num_trees + 3)), "Computing edge weights.")
    weights = edge_weights(g, hub_penalty)
    
    # Find first steiner trees
    seeds = list(filter(lambda s: s in seed_map, seeds))
    task_hook.set_progress(2 / (float(num_trees + 3)), "Computing Steiner tree 1 of {}.".format(num_trees))
    first_tree = steiner_tree(g, seeds, seed_map, weights, hub_penalty > 0)
    num_found_trees = 1
    tree_edges = []
    for tree_edge in first_tree.edges():
        source_name = first_tree.vertex_properties[node_name_attribute][first_tree.vertex_index[tree_edge.source()]]
        target_name = first_tree.vertex_properties[node_name_attribute][first_tree.vertex_index[tree_edge.target()]]
        tree_edges.append((gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute], match=source_name)[0], gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute], match=target_name)[0]))
    cost_first_tree = sum([weights[g.edge(source, target)] for source, target in tree_edges])
    returned_nodes = set(int(gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute], match=first_tree.vertex_properties[node_name_attribute][node])[0]) for node in range(first_tree.num_vertices()))
    if num_trees > 1:
        is_bridge = find_bridges(g)
        edge_filter = g.new_edge_property("boolean", True)
        found_new_tree = True
        while len(tree_edges) > 0:
            if found_new_tree:
                task_hook.set_progress(float(num_found_trees + 2) / (float(num_trees + 3)), "Computing Steiner tree {} of {}.".format(num_found_trees + 1, num_trees))
            found_new_tree = False
            tree_edge = tree_edges.pop()
            g_edge = g.edge(tree_edge[0], tree_edge[1])
            if not is_bridge[g_edge]:
                edge_filter[g_edge] = False
                g.set_edge_filter(edge_filter)
                next_tree = steiner_tree(g, seeds, seed_map, weights, hub_penalty > 0)
                next_tree_edges = set()
                for next_tree_edge in next_tree.edges():
                    source_name = next_tree.vertex_properties[node_name_attribute][next_tree.vertex_index[next_tree_edge.source()]]
                    target_name = next_tree.vertex_properties[node_name_attribute][next_tree.vertex_index[next_tree_edge.target()]]
                    next_tree_edges.add((gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute], match=source_name)[0],gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute], match=target_name)[0]))
                cost_next_tree = sum([weights[g.edge(source, target)] for source, target in next_tree_edges])
                if cost_next_tree <= cost_first_tree * ((100.0 + tolerance) / 100.0):
                    found_new_tree = True
                    num_found_trees += 1
                    for node in range(next_tree.num_vertices()):
                        returned_nodes.add(int(gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute],match=next_tree.vertex_properties[node_name_attribute][node])[0]))
                    removed_edges = []
                    for source, target in tree_edges:
                        if not ((source, target) in set(next_tree_edges)) or ((target, source) in set(next_tree_edges)):
                            removed_edges.append((source, target))
                    for edge in removed_edges:
                        tree_edges.remove(edge)
                g.clear_filters()
                edge_filter[g_edge] = True
            if num_found_trees >= num_trees:
                break
    task_hook.set_progress((float(num_trees + 2)) / (float(num_trees + 3)), "Formatting results")
    returned_edges = []
    for node in returned_nodes:
        for neighbor in g.get_all_neighbors(node):
            if int(neighbor) > node and int(neighbor) in returned_nodes:
                returned_edges.append((node, int(neighbor)))

    accepted_nodes = [g.vertex_properties[node_name_attribute][node] for node in returned_nodes]
    accepted_nodes_without_seeds = [g.vertex_properties[node_name_attribute][node] for node in returned_nodes if node not in seed_ids]
    
    # avoid duplicate edges and sort edge source/target lexicographically
    edges_unique = set()
    for source, target in returned_edges:
        a = g.vertex_properties[node_name_attribute][source]
        b = g.vertex_properties[node_name_attribute][target]
        if a > b:
            # flip source and target to get a < b
            tmp = a
            a = b
            b = tmp
        edges_unique.add((a, b))
    subgraph = {"nodes": accepted_nodes,
                "edges": [{"from": source, "to":target} for
                          source, target in edges_unique]}
    
    node_types = {g.vertex_properties[node_name_attribute][node]: g.vertex_properties["type"][node] for node in returned_nodes}
    is_seed = {g.vertex_properties[node_name_attribute][node]: node in set(seed_ids) for node in returned_nodes}
    
    task_hook.set_results({
        "network": subgraph,
        "node_attributes": {"node_types": node_types, "is_seed": is_seed},
        "target_nodes": accepted_nodes_without_seeds,
        'gene_interaction_dataset': ppi_dataset,
        'drug_interaction_dataset': pdi_dataset
    })
