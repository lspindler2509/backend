from drugstone.util.property_calulations import calculate_properties
from drugstone.util.query_db import query_proteins_by_identifier, fetch_edges_for_proteins
from tasks.util.custom_network import add_edges, remove_ppi_edges
from tasks.task_hook import TaskHook
import graph_tool as gt
from drugstone.models import *
from drugstone.serializers import *
import os


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
    
    no_default_edges = task_hook.parameters.get("exclude_drugstone_ppi_edges", False)
    
    
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
    if task_hook.parameters["config"].get("reviewed", False):
        filename += "_reviewed"
    filename = os.path.join(task_hook.data_directory, filename + ".gt")
    g = gt.load_graph(filename)
    if custom_edges:
        if no_default_edges:
          # clear all edges with type "protein-protein"
          g = remove_ppi_edges(g)
        edges = task_hook.parameters.get("input_network")['edges']
        g = add_edges(g, edges)
        
    task_hook.set_progress(2 / 4.0, "Get all first neighbors from database.")
    
    # Map seeds to UniProt IDs
    reviewed = task_hook.parameters["config"].get("reviewed", False)
    seeds_mapped, _ = query_proteins_by_identifier(set(seeds), identifier_key, reviewed)
    seeds_uniprot = {uniprot for node in seeds_mapped if node.get("uniprot") for uniprot in node["uniprot"]}
    
    
    # Fetch all edges where seeds are either from_protein or to_protein
    interaction_objects = fetch_edges_for_proteins(
        ppi_dataset['name'], 
        ppi_dataset['licenced'], 
        seeds_uniprot
    )
    
    # Extract all unique proteins from these edges (seeds + first neighbors)
    # Start with seeds_uniprot to ensure seeds are always included
    all_protein_ids = set(seeds_uniprot)
    for interaction in interaction_objects:
        all_protein_ids.add(interaction.from_protein.uniprot_code)
        all_protein_ids.add(interaction.to_protein.uniprot_code)
    
    # Map all proteins (seeds + first neighbors) - use "uniprot" since all_protein_ids contains UniProt codes
    all_nodes_mapped, _ = query_proteins_by_identifier(all_protein_ids, "uniprot", reviewed)
    
    
    # Create mapping from identifier_key to node
    nodes_mapped_dict = {}
    for node in all_nodes_mapped:
        if node.get(identifier_key):
            nodes_mapped_dict[node[identifier_key][0]] = node
    
    # Create drugstone_mapping: drugstone_id -> identifier_key value
    drugstone_mapping = {}
    for node in all_nodes_mapped:
        if node.get("drugstone_id") and node.get(identifier_key):
            drugstone_mapping[node["drugstone_id"][0]] = node[identifier_key][0]
    
    # Create node details
    all_nodes_mapped_list = []
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
        isReviewed = nodes_mapped_dict[node].get("is_reviewed", False)
        if node in set(seedSet):
            isSeed[node] = True
            group = "seedNode"
        else:
            isSeed[node] = False
            group = "firstNeighbor"
            
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
            "groupId": group,
            "cellular_component": cellular_component,
            "layer": layer,
            "isReviewed": isReviewed,
        }
        all_nodes_mapped_list.append(mapped_node)
    
    task_hook.set_progress(3 / 4.0, "Get all edges within the subnetwork.")
    
    # Get all UniProt codes in the subnetwork (seeds + first neighbors)
    # Include all UniProt IDs, not just the first one
    subnet_uniprot = set(seeds_uniprot)  # Start with seeds to ensure they're always included
    for node in all_nodes_mapped:
        if node.get("uniprot"):
            subnet_uniprot.update(node["uniprot"])  # Add all UniProt IDs
    
    # Fetch ALL edges between subnetwork nodes (not just those involving seeds)
    # This includes edges between first neighbors that don't involve seeds directly
    # Use require_both_nodes=True for more efficient query
    subnet_interactions = fetch_edges_for_proteins(
        ppi_dataset['name'], 
        ppi_dataset['licenced'], 
        subnet_uniprot,
        require_both_nodes=True
    )
    
    # Convert edges directly from interaction objects (not serialized)
    # This way we have access to the actual Protein objects with uniprot_code
    edges = []
    uniprot_to_identifier = {}
    # Create mapping from all UniProt IDs to identifier_key value
    for node in all_nodes_mapped:
        if node.get("uniprot") and node.get(identifier_key):
            identifier_value = node[identifier_key][0]
            for uniprot in node["uniprot"]:
                uniprot_to_identifier[uniprot] = identifier_value
    
    for interaction in subnet_interactions:
        from_uniprot = interaction.from_protein.uniprot_code
        to_uniprot = interaction.to_protein.uniprot_code
        
        if from_uniprot in uniprot_to_identifier and to_uniprot in uniprot_to_identifier:
            edge = {
                "from": uniprot_to_identifier[from_uniprot],
                "to": uniprot_to_identifier[to_uniprot],
                "is_directed": interaction.is_directed,
                "is_stimulation": interaction.is_stimulation,
                "is_inhibition": interaction.is_inhibition,
                "dataset": ppi_dataset['name'],
            }
            edges.append(edge)

    # Filter nodes to keep only upstream regulators if parameter is set
    only_upstream_regulators = task_hook.parameters.get("only_upstream_regulators", False)
    is_omnipath = ppi_dataset["name"] == "OmniPath"
    if only_upstream_regulators:
        # Only allow upstream regulator filtering for OmniPath (which has directed edges)
        if not is_omnipath:
            raise ValueError("Filtering for upstream regulators is only supported for OmniPath dataset (which has directed edges).")
        node_ids = {node["id"] for node in all_nodes_mapped_list}
        seed_ids = {node["id"] for node in all_nodes_mapped_list if node.get("group") == "seedNode"}
        
        # Track nodes with outgoing directed edges to seeds
        # Simple logic: if a node has a directed edge to ANY seed (where seed is target), keep it
        upstream_regulators = set()
        
        for edge in edges:
            from_node = edge.get("from")
            to_node = edge.get("to")
            is_directed = edge.get("is_directed", False)
            
            if isinstance(is_directed, str):
                is_directed = is_directed.lower() == "true"
            
            # Simple check: is it a directed edge where target is a seed and source is not a seed?
            if is_directed and to_node in seed_ids and from_node in node_ids and from_node not in seed_ids:
                upstream_regulators.add(from_node)
        
        # Determine which nodes to keep:
        # - Seeds always stay
        # - Nodes with outgoing directed edges to seeds (upstream regulators) stay
        nodes_to_keep = seed_ids | upstream_regulators
        
        # Filter nodes and edges BEFORE calculating properties
        all_nodes_mapped_list = [node for node in all_nodes_mapped_list if node["id"] in nodes_to_keep]
        edges = [edge for edge in edges 
                 if edge.get("from") in nodes_to_keep and edge.get("to") in nodes_to_keep]
        
        # Update nodes_mapped_dict after filtering
        nodes_mapped_dict = {}
        for node in all_nodes_mapped_list:
            if node.get(identifier_key):
                nodes_mapped_dict[node[identifier_key][0]] = node
    
    edges_for_properties = []
    for edge in edges:
        from_id = edge["from"]
        to_id = edge["to"]
        if from_id in nodes_mapped_dict and to_id in nodes_mapped_dict:
            edges_for_properties.append({
                "from": from_id,
                "to": to_id
            })
    
    all_nodes_mapped = calculate_properties(all_nodes_mapped_list, g, identifier_key, edges_for_properties, True)

    # Automatic SPD cutoff suggestion if network exceeds 250 nodes
    MAX_NODES = 250
    network_initial = {"nodes": all_nodes_mapped, "edges": edges}
    
    if len(all_nodes_mapped) > MAX_NODES:
        def get_spd_value(node):
            spd = node.get("properties", {}).get("spd")
            if spd is None:
                return float('inf')
            try:
                return float(spd)
            except (ValueError, TypeError):
                return float('inf')
        
        # Group ALL nodes by SPD value (including seeds, which have SPD = 1)
        from collections import defaultdict
        nodes_by_spd = defaultdict(list)
        for node in all_nodes_mapped:
            spd = get_spd_value(node)
            nodes_by_spd[spd].append(node)
        print(nodes_by_spd, "\n")
        
        # Sort SPD values descending (higher = closer to seeds, seeds have SPD = 1)
        sorted_spd_values = sorted([spd for spd in nodes_by_spd.keys() if spd != float('inf')], reverse=True)
        print(sorted_spd_values, "\n")
        
        nodes_to_keep_ids = set()
        nodes_to_keep_count = 0
        suggested_cutoff = None
        
        for spd in sorted_spd_values:
            count_with_this_spd = len(nodes_by_spd[spd])
            if nodes_to_keep_count + count_with_this_spd < MAX_NODES:
                for node in nodes_by_spd[spd]:
                    nodes_to_keep_ids.add(node["id"])
                nodes_to_keep_count += count_with_this_spd
                suggested_cutoff = spd
            else:
                break
        
        if suggested_cutoff is not None:
            print(suggested_cutoff, "\n")
            all_nodes_mapped = [node for node in all_nodes_mapped if node["id"] in nodes_to_keep_ids]
            edges = [edge for edge in edges 
                     if edge.get("from") in nodes_to_keep_ids and edge.get("to") in nodes_to_keep_ids]
            task_hook.parameters["suggested_spd_cutoff"] = suggested_cutoff
        else:
            print("No suggested cutoff found")

    task_hook.set_progress(4 / 4.0, "Returning results.")

    result = {"nodes": all_nodes_mapped, "edges": edges}
    task_hook.parameters["algorithm"] = "first-neighbor"
    
    # Store cutoff at result level (same as when user manually prunes)
    result_dict = {
        "algorithm": "first_neighbor",
        "network": result,
        "network_initial": network_initial,
        "parameters": task_hook.parameters,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
        "node_attributes": {"is_seed": isSeed},
    }
    
    # If automatic cutoff was applied, store it at result level (like manual pruning)
    if task_hook.parameters.get("suggested_spd_cutoff") is not None:
        result_dict["cutoff"] = task_hook.parameters["suggested_spd_cutoff"]
        result_dict["pruneOrphanNodes"] = False
        result_dict["automaticCutoff"] = True
    
    task_hook.set_results(result_dict)
