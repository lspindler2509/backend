from drugstone.util.query_db import query_proteins_by_identifier
from tasks.util.custom_network import add_edges
from tasks.task_hook import TaskHook
import graph_tool as gt
from drugstone.models import *
from drugstone.serializers import *
import os
import graph_tool.util as gtu



def first_neighbor(task_hook: TaskHook):
    r"""
    Get all first neighbors of the seed genes.

    Parameters
    ----------
    

    num_threads : int, optional (default: 1)
      Number of threads. Requires that graph_tool is compiled with OpenMP support.
      Should not be exposed in the frontend.
      
    Returns
    -------
    results : {
        "algorithm": "first_neighbor", # Name of the algorithm.
        "network":result, # The network with the seed nodes and their first neighbors.
        "parameters": task_hook.parameters,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
        "node_attributes":
            {
                "is_seed": isSeed,
            },
    }

    Notes
    -----
    This implementation is based on graph-tool, a very efficient Python package for network
    analysis with C++ backend and multi-threading support. Installation instructions for graph-tool
    can be found at https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions.

    """
    

    # Type: list of str
    # Semantics: Names of the seed proteins. Use UNIPROT IDs for host proteins, and
    #            names of the for SARS_CoV2_<IDENTIFIER> (e.g., SARS_CoV2_ORF6) for
    #            virus proteins.
    # Reasonable default: None, has to be selected by user via "select for analysis"
    #            utility in frontend.
    # Acceptable values: UNIPROT IDs, identifiers of viral proteins.
    seeds = task_hook.parameters["seeds"]


    # Type: int.
    # Semantics: Number of threads used for running the analysis.
    # Example: 1.
    # Reasonable default: 1.
    # Note: We probably do not want to expose this parameter to the user.
    num_threads = task_hook.parameters.get("num_threads", 1)

    ppi_dataset = task_hook.parameters.get("ppi_dataset")

    pdi_dataset = task_hook.parameters.get("pdi_dataset")

    id_space = task_hook.parameters["config"].get("identifier", "symbol")

    custom_edges = task_hook.parameters.get("custom_edges", False)
    
    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    identifier_key = id_space
    if id_space == "ncbi":
        identifier_key = "entrez"
    elif id_space == "ensembl":
        identifier_key = "ensg"
    
    # Parsing input file.
    task_hook.set_progress(1 / 4.0, "Parsing input.")
    
    filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset['licenced'] or pdi_dataset['licenced']:
        filename += "_licenced"
    filename = os.path.join(task_hook.data_directory, filename + ".gt")
    g = gt.load_graph(filename)
    if custom_edges:
        edges = task_hook.parameters.get("input_network")['edges']
        g = add_edges(g, edges)
        
    task_hook.set_progress(2 / 4.0, "Get all first neighbors.")
    
    node_name_attribute = "internal_id"
    node_mapping = {}
    node_mapping_reverse = {}
    all_neighbors_numbers = set()
    for seed in seeds:
        found = gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute], match=seed)
        if len(found) > 0:
            found_node = int(found[0])
            node_mapping[seed] = found_node
            node_mapping_reverse[found_node] = seed
            new_neighbors = g.get_all_neighbors(found_node)
            for neighbor in new_neighbors:
                # We only want to add the neighbor if it is a protein.
                if g.vertex_properties['drug_id'][neighbor] == "":
                    all_neighbors_numbers.add(neighbor)
                    node = g.vertex_properties[node_name_attribute][neighbor]
                    node_mapping[node] = neighbor
                    node_mapping_reverse[neighbor] = node
    
    all_neighbors = list(node_mapping.keys())
    
    nodes_mapped, identifier = query_proteins_by_identifier(all_neighbors, identifier_key)
    nodes_mapped_dict = {node[identifier][0]: node for node in nodes_mapped}
    
    # Get the node details.
    all_nodes_mapped = []
    isSeed = {}
    seedSet = set(seeds)
    for node in nodes_mapped_dict.keys():
        drugstone_id = nodes_mapped_dict[node]["drugstone_id"]
        uniprot = nodes_mapped_dict[node]["uniprot"]
        symbol = nodes_mapped_dict[node].get("symbol", "")
        protein_name = nodes_mapped_dict[node]["protein_name"]
        entrez = nodes_mapped_dict[node]["entrez"]
        cellular_component = nodes_mapped_dict[node].get("cellular_component", [])
        layer = nodes_mapped_dict[node].get("layer", "")
        ensg = nodes_mapped_dict[node].get("ensg", "")
        if node in set(seedSet):
            isSeed[node] = True
            group = "seedNode"
        else:
            isSeed[node] = False
            group = "default"
            
        mapped_node = {
            "id": nodes_mapped_dict[node][identifier_key][0],
            "drugstone_id": drugstone_id,
            "drugstone_type": "protein",
            "uniprot": uniprot,
            "symbol": symbol,
            "protein_name": protein_name,
            "entrez": entrez,
            "ensg": ensg,
            "label": nodes_mapped_dict[node][identifier_key][0],
            "group": group,
            "cellular_component": cellular_component,
            "layer": layer,
        }
        all_nodes_mapped.append(mapped_node) 
        
    task_hook.set_progress(3 / 4.0, "Get all edges connecting the seed nodes with the first neighbors.")
            
    edges_unique = set()
    for node in node_mapping.keys():
        for neighbor in g.get_all_neighbors(node_mapping[node]):
            if int(neighbor) > int(node_mapping[node]) and int(neighbor) in all_neighbors_numbers:
                first_key = next(iter(node_mapping_reverse))
                if isinstance(first_key, int):
                    neighbor_key = int(neighbor)
                else:
                    neighbor_key = str(int(neighbor))
                edges_unique.add((node, node_mapping_reverse[neighbor_key]))
    
    edges = [{"from": source, "to":target} for source, target in edges_unique]


    # return the results.
    task_hook.set_progress(4 / 4.0, "Returning results.")
    
    result = {
        "nodes": all_nodes_mapped,
        "edges": edges,
    }
    task_hook.parameters["algorithm"] = "first-neighbor"
    task_hook.set_results({
        "algorithm": "first_neighbor",
        "network":result,
        "parameters": task_hook.parameters,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
        "node_attributes":
            {
                "is_seed": isSeed,
            },
    })
