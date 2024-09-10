from tasks.util.custom_network import add_edges
from tasks.task_hook import TaskHook
import graph_tool as gt
from drugstone.models import *
from drugstone.serializers import *
import os
import networkx as nx
import community as community_louvain
from collections import Counter
import matplotlib.colors as mcolors
import graph_tool.util as gtu
import leidenalg
from igraph import Graph

import distinctipy

def generate_color_palette(num_colors):
    palette = distinctipy.get_colors(num_colors)

    hex_colors = [mcolors.rgb2hex(color) for color in palette]

    return hex_colors


def add_cluster_groups_to_config(config, partition):
    unique_clusters = set(partition.values())
    color_palette = generate_color_palette(len(unique_clusters))
    
    for cluster_id, color in zip(unique_clusters, color_palette):
        group_name = f"Cluster {cluster_id}"
        group_id = f"cluster{cluster_id}"
        
        if not config["node_groups"].get(group_id):
            config["node_groups"][group_id] = {
                "group_name": group_name,
                "color": {
                    "border": color,
                    "background": color,
                    "highlight": {
                        "border": color,
                        "background": color
                    }
                },
                "shape": "circle",
                "type": "gene",
                "border_width": 0,
                "border_width_selected": 0,
                "font": {
                    "color": "#FFFFFF" if is_dark_color(color) else "#000000",  # Setze Schriftfarbe auf Weiß, wenn Hintergrund dunkel ist
                    "size": 14,
                    "face": "arial",
                    "stroke_width": 0,
                    "stroke_color": "#ffffff",
                    "align": "center",
                    "bold": False,
                    "ital": False,
                    "boldital": False,
                    "mono": False
                },
                "shadow": True,
                "group_id": group_id
            }
    return config

def is_dark_color(color):
    r, g, b, _ = mcolors.to_rgba(color)
    brightness = (r * 299 + g * 587 + b * 114) / 1000
    return brightness < 0.5


def leiden_clustering(task_hook: TaskHook):
    r"""
    Perform Leiden clustering.

    Parameters
    ----------
    

    num_threads : int, optional (default: 1)
      Number of threads. Requires that graph_tool is compiled with OpenMP support.
      Should not be exposed in the frontend.
      
    ignore_isolated : bool, optional (default: True)
    Specifys to not include the isolated nodes in the clustering results, since they build their own clusteranyways.
    They will keep their old group, so no distinct colours are wasted.


    Returns
    -------
     results : {
        "algorithm": "louvain_clustering", # Name of the algorithm.
        "network":result, # The network with the clustering results, the group specifies the cluster.
        "table_view": table_view_results, # some statistics about the clustering results.
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

    References
    ----------
    ..  [1] Traag, Vincent, Fabio Zanini, Ryan Gibson, Oren Ben-Kiki, Tom Kelly, Britny Farahdel, and Dafne van Kuppevelt. “Vtraag/Leidenalg: 0.10.0.” Zenodo, July 14, 2023. https://doi.org/10.5281/zenodo.8147844.
        [2] Antonov, Michael, Gábor Csárdi, Szabolcs Horvát, Kirill Müller, Tamás Nepusz, Daniel Noom, Maëlle Salmon, Vincent Traag, Brooke Foucault Welles, and Fabio Zanini. “Igraph Enables Fast and Robust Network Analysis across Programming Languages.” arXiv, November 16, 2023. https://doi.org/10.48550/arXiv.2311.10260.
        [3] Roberts, Jack, Jon Crall, Kian-Meng Ang, and Yannick Brandt. “Alan-Turing-Institute/Distinctipy: V1.3.4.” Zenodo, January 10, 2024. https://doi.org/10.5281/zenodo.10480933.

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
    
    ignore_isolated = task_hook.parameters.get("ignore_isolated", True)
    
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
        
    node_name_attribute = "internal_id"
    node_mapping = {}
    node_mapping_reverse = {}
    for seed in seeds:
        found = gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute], match=seed)
        if len(found) > 0:
            found_node = int(found[0])
            node_mapping[seed] = found_node
            node_mapping_reverse[found_node] = seed
        
    all_nodes_int = set([int(node_mapping[gene]) for gene in seeds if gene in node_mapping])
    edges_unique = set()
    for node in node_mapping.keys():
        for neighbor in g.get_all_neighbors(node_mapping[node]):
            if int(neighbor) > int(node_mapping[node]) and int(neighbor) in all_nodes_int:
                first_key = next(iter(node_mapping_reverse))
                if isinstance(first_key, int):
                    neighbor_key = int(neighbor)
                else:
                    neighbor_key = str(int(neighbor))
                edges_unique.add((node, node_mapping_reverse[neighbor_key]))

    
    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    edges = [{"from": source, "to":target} for source, target in edges_unique]
    nodes = task_hook.parameters.get("input_network")['nodes']
        
    isSeed = {}
    seedSet = set(seeds)
    g = Graph(directed=False)
    for node in nodes:
        if node["id"] in seedSet:
            isSeed[node["id"]] = True
            g.add_vertex(name=node["id"])
    for edge in edges:
        if edge["from"] in seedSet and edge["to"] in seedSet:
            g.add_edge(edge["from"], edge["to"])
            
    if ignore_isolated:
        isolated_nodes = [v.index for v in g.vs if v.degree() == 0]
        g.delete_vertices(isolated_nodes)


    task_hook.set_progress(3 / 4.0, "Perform Louvain clustering.")

    partition: leidenalg.VertexPartition.ModularityVertexPartition = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)
    
    partition_dict = {}
    counter = 0
    for cluster in partition._formatted_cluster_iterator():
        nodes_cluster = cluster.split(", ")
        for node in nodes_cluster:
            partition_dict[node] = counter
        counter += 1
                        
    task_hook.set_progress(3 / 4.0, "Parse clustering results.")

    config = add_cluster_groups_to_config(task_hook.parameters["config"], partition_dict)
    task_hook.parameters["config"] = config
   
    table_view_results = []
    total_nodes = len(partition_dict)
    cluster_counts = Counter(partition_dict.values())
    for cluster_id, count in cluster_counts.items():
        percentage = round((count / total_nodes) * 100, 2)
        table_view_results.append({
            "cluster_id": cluster_id,
            "count": str(count) + "/"+ str(total_nodes),
            "percentage": percentage,
        })
    
    filtered_nodes = []
    for node in nodes:
        if node["id"] in seedSet:
            if node["id"] in partition_dict:
                cluster = partition_dict[node["id"]]
                group_id = f"cluster{cluster}"
                node["group"] = group_id
                node["cluster"] = str(cluster)
                filtered_nodes.append(node)
            else:
                # node in seeds but was isolated and those are ignored -> keep old group
                node["cluster"] = "none"
                filtered_nodes.append(node)
                


    # return the results.
    task_hook.set_progress(4 / 4.0, "Returning results.")
    
    result = {
        "nodes": filtered_nodes,
        "edges": edges,
    }
    task_hook.parameters["algorithm"] = "leiden-clustering"
    task_hook.set_results({
        "algorithm": "leiden_clustering",
        "network":result,
        "table_view": table_view_results,
        "parameters": task_hook.parameters,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
        "node_attributes":
            {
                "is_seed": isSeed,
            },
    })
