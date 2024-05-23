import numpy as np
from tasks.util.custom_network import add_edges, remove_ppi_edges, filter_proteins
from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.scores_to_results import scores_to_results
from tasks.util.edge_weights import edge_weights
from tasks.task_hook import TaskHook
import graph_tool as gt
import graph_tool.topology as gtt
import os.path
import sys

def closeness_centrality(task_hook: TaskHook):
    r"""Computes closeness centrality w.r.t. seed nodes.
    
    The closeness centrality of a node :math:`u` in a graph :math:`G=(V,E)` is defined 
    as :math:`|V|/\sum_{v\in V}d(v,u)`, where :math:`d(v,u)` is the length of a shortest path 
    from :math:`v` to :math:`u` [1_]. We here use the modified version :math:`|S|/\sum_{v\in S}d(v,u)`
    suggested in [2_], where :math:`S\subseteq V` is a set of selected seed nodes.
    
    Parameters
    ----------
    # strain_or_drugs : str
    #   The virus strain for which the analysis should be run, or the string literal "drugs"
    #   (if used for ranking drugs).
    target_or_drugs : str
      The string literal "target", or the string literal "drugs"
      (if used for ranking drugs).
    seeds : list of str
      A list of UNIPROT IDs identifying the seed proteins.
      
    datasets : list of str
      List of datasets whose nodes and returned_edges should be considered for the analysis.
      If empty, all available datasets are employed.
      
    ignored_edge_types : list of str, optional (default: [])
      Edges whose types are contained in this list are ignored in the analysis.
      
    include_indirect_drugs : bool, optional (default: False)
      If True, also drugs targeting interactors of the seeds are considered when ranking drugs.
      
    include_non_approved_drugs : bool, optional (default: False)
      If True, also non-approved drugs are considered when ranking drugs.
      
    result_size : int, optional (default: 20)
      The number of new candidate proteins to be returned.
      
    num_threads : int, optional (default: 1)
      Number of threads. Requires that graph_tool is compiled with OpenMP support.
      Should not be exposed in the frontend.
      
    Returns
    -------
    results : {"networks": list of dict, "node_attributes": list of dict}
      "networks": A one-element list containing the subgraph induced by the result 
        proteins and the seeds, along with the maximal score of all nodes, the maximal
        score off all nodes except virus proteins, and the maximal scores of all nodes
        except seed nodes.
      "node_attributes": A one-element list containing a dictionary with the following 
        attributes for all nodes in the returned network:
        "node_types": The type of the nodes (either "virus", "host", or "drug").
        "is_seed": A flag that specifies whether the node is a seed.
        "scores": The un-normalized scores for all non-seed nodes (nan for the virus proteins).
      
    Notes
    -----
    This implementation is based on graph-tool, a very efficient Python package for network
    analysis with C++ backend and multi-threading support. Installation instructions for graph-tool
    can be found at https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions.
    
    References
    ----------
    .. [1] G. Sabidussi, The Centrality Index of a Graph, Psychometrika 31(4), 1966, pp. 581–603,
       https://doi.org/10.1007/bf02289527.
       [2] T. Kacprowski, N.T. Doncheva, M. Albrecht, NetworkPrioritizer: A Versatile Tool for 
       Network-Based Prioritization of Candidate Disease Genes or Other Molecules, Bioinformatics 29(11),
       pp. 1471-1473, https://doi.org/10.1093/bioinformatics/btt164.  
    """
    
    # Type: list of str
    # Semantics: Names of the seed proteins. Use UNIPROT IDs for host proteins, and 
    #            names of the for SARS_CoV2_<IDENTIFIER> (e.g., SARS_CoV2_ORF6) for
    #            virus proteins.
    # Reasonable default: None, has to be selected by user via "select for analysis"
    #            utility in frontend.
    # Acceptable values: UNIPROT IDs, identifiers of viral proteins.
    seeds = task_hook.parameters["seeds"]
    
    # Type: str.
    # Semantics: The virus strain for which the analysis should be run, or the 
    #            string literal "drugs" (if used for ranking drugs).
    # Example: "SARS_CoV2"
    # Reasonable default: None, has to be specified by the caller.
    ### Acceptable values: "SARS_CoV2", "drugs", ...
    # strain_or_drugs = task_hook.parameters.get("strain_or_drugs", "SARS_CoV2")
    # Acceptable values: "PPI", "PPDr"
    # target_or_drugs = task_hook.parameters.get("target_or_drugs", "PPI")
    
    # Type: list of str.
    # Semantics: The datasets which should be considered for the analysis.
    # Example: ["Krogan", "TUM"].
    # Note: If empty, all available datasets are used. When ranking drugs, the 
    #       default [] should be used.
    # Reasonable default: [].
    # Acceptable values: "Krogan", "TUM".
    # datasets = task_hook.parameters.get("datasets", [])
    
    # Type: list of str.
    # Semantics: Virus-host edge types which should be ignored for the analysis.
    # Example: ["Overexpression"].
    # Note: If empty, all available edge types are used. When ranking drugs, the 
    #       default [] should be used.
    # Reasonable default: [].
    # Acceptable values: "AP-MS", "overexpression".
    # ignored_edge_types = task_hook.parameters.get("ignored_edge_types", [])
    
    # Type: bool
    # Semantics: Sepcifies whether also drugs targeting interactors of the seeds should be considered.
    # Example: False.
    # Reasonable default: False.
    # Has no effect unless trust_rank.py is used for ranking drugs.
    include_indirect_drugs = task_hook.parameters.get("include_indirect_drugs", False)
    # Type: bool
    # Semantics: Sepcifies whether should be included in the analysis when ranking drugs.
    # Example: False.
    # Reasonable default: False.
    # Has no effect unless trust_rank.py is used for ranking drugs.
    include_non_approved_drugs = task_hook.parameters.get("include_non_approved_drugs", False)
    
    # Type bool.
    # Semantics: Ignore viral proteins which are not selected as seeds.
    # Example: False.
    # Reasonable default: False.
    # Has no effect when the algorithm is used for ranking drugs.
    # ignore_non_seed_baits = task_hook.parameters.get("ignore_non_seed_baits", False)
    
    # Type: int.
    # Semantics: Number of returned proteins.
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

    search_target = task_hook.parameters.get("target", "drug")

    filterPaths = task_hook.parameters.get("filter_paths", True)
    
    # Parsing input file.
    task_hook.set_progress(0 / 4.0, "Parsing input.")

    id_space = task_hook.parameters["config"].get("identifier", "symbol")
    
    custom_edges = task_hook.parameters.get("custom_edges", False)
    
    no_default_edges =    no_default_edges = task_hook.parameters.get("exclude_drugstone_ppi_edges", False)
    
    custom_nodes = task_hook.parameters.get("network_nodes", False)

    filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset['licenced'] or pdi_dataset['licenced']:
        filename += "_licenced"
    filename = os.path.join(task_hook.data_directory, filename + ".gt")
    g, seed_ids, drug_ids = read_graph_tool_graph(filename, seeds, id_space, max_deg, include_indirect_drugs, include_non_approved_drugs, search_target)
    
    if custom_edges:
      if no_default_edges:
        # clear all edges with type "protein-protein"
        g = remove_ppi_edges(g)
      g = add_edges(g, custom_edges)
    
    if custom_nodes:
      # remove all nodes with internal_id not in custom_nodes from g
      g, seed_ids, drug_ids = filter_proteins(g, custom_nodes, drug_ids, seeds)
    
    task_hook.set_progress(1 / 4.0, "Computing edge weights.")
    weights = edge_weights(g, hub_penalty) 
    
    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    # Call graph-tool to compute TrustRank.
    task_hook.set_progress(2 / 4.0, "Computing shortest path closeness centralities.")
    all_dists = []
    # score_nodes = drug_ids if search_target == 'drug' else seed_ids
    for node in seed_ids:
        dists = gtt.shortest_distance(g, node, weights=weights).get_array()
        dists[dists == np.inf] = 99999999999
        all_dists.append(dists+1)

    scores = len(seed_ids) / (sum([dists for dists in all_dists]))

    # Compute and return the results.
    task_hook.set_progress(3 / 4.0, "Formatting results.")
    
    task_hook.set_results(scores_to_results(search_target, result_size, g, seed_ids, drug_ids, scores, ppi_dataset, pdi_dataset, filterPaths))

