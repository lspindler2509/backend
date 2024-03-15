from tasks.util.custom_edges import add_edges
from tasks.util.read_graph_tool_graph import read_graph_tool_graph
from tasks.util.scores_to_results import scores_to_results
from tasks.util.edge_weights import edge_weights
from tasks.task_hook import TaskHook
import graph_tool as gt
import graph_tool.topology as gtt
import sys
import gseapy as gp
from drugstone.models import *
from drugstone.serializers import *


def pathway_enrichment(task_hook: TaskHook):
    r"""
    Perform pathway enrichment analysis.

    Parameters
    ----------

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
    .. [1] L. Freeman, A Set of Measures of Centrality Based on Betweenness, Sociometry 40(1), 
       1977, pp. 35–41, https://doi.org/10.2307/3033543.
       [2] T. Kacprowski, N.T. Doncheva, M. Albrecht, NetworkPrioritizer: A Versatile Tool for 
       Network-Based Prioritization of Candidate Disease Genes or Other Molecules, Bioinformatics 29(11),
       2013, pp. 1471-1473, https://doi.org/10.1093/bioinformatics/btt164.  
    """

    # Type: list of str
    # Semantics: Names of the seed proteins. Use UNIPROT IDs for host proteins, and
    #            names of the for SARS_CoV2_<IDENTIFIER> (e.g., SARS_CoV2_ORF6) for
    #            virus proteins.
    # Reasonable default: None, has to be selected by user via "select for analysis"
    #            utility in frontend.
    # Acceptable values: UNIPROT IDs, identifiers of viral proteins.
    seeds = task_hook.parameters["seeds"]

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

    # Type: int.
    # Semantics: Number of returned proteins.
    # Example: 20.
    # Reasonable default: 20.
    # Acceptable values: integers n with n > 0.
    result_size = task_hook.parameters.get("result_size", 20)

    # Type: int.
    # Semantics: Number of threads used for running the analysis.
    # Example: 1.
    # Reasonable default: 1.
    # Note: We probably do not want to expose this parameter to the user.
    num_threads = task_hook.parameters.get("num_threads", 1)

    # ppi_dataset = task_hook.parameters.get("ppi_dataset")

    # pdi_dataset = task_hook.parameters.get("pdi_dataset")

    search_target = task_hook.parameters.get("target", "drug-target")

    filterPaths = task_hook.parameters.get("filter_paths", True)

    id_space = task_hook.parameters["config"].get("identifier", "symbol")

    custom_edges = task_hook.parameters.get("custom_edges", False)

    # Parsing input file.
    task_hook.set_progress(0 / 3.0, "Parsing input.")
    # filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    # if ppi_dataset['licenced'] or pdi_dataset['licenced']:
    #     filename += "_licenced"
    # filename = os.path.join(task_hook.data_directory, filename + ".gt")
    # g, seed_ids, drug_ids = read_graph_tool_graph(
    #     filename,
    #     seeds,
    #     id_space,
    #     max_deg,
    #     include_indirect_drugs,
    #     include_non_approved_drugs,
    #     target=search_target
    # )

    # if custom_edges:
    #     edges = task_hook.parameters.get("input_network")['edges']
    #     g = add_edges(g, edges)

    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)

    # Calculate pathway enrichment
    if id_space == "symbol":
        seeds_symbol = seeds
    else:
        seeds_symbol = []
        for seed in seeds:
            if id_space == "uniprot":
                seed_without_identifier = seed.replace("uniprot.", "")
                protein = models.Protein.objects.filter(
                    uniprot_code=seed_without_identifier
                ).last()
                seeds_symbol.append(protein.gene)
            if id_space == "entrez":
                protein = models.Protein.objects.filter(
                    entrez=seed
                ).last()
                seeds_symbol.append(protein.gene)

    print(seeds_symbol)

    enr = gp.enrichr(gene_list=seeds_symbol,
                     gene_sets=['WikiPathway_2023_Human', 'KEGG_2021_Human', 'Reactome_2022'],
                     organism='human',
                     outdir=None,
                     )
    print(enr.results.head(10))

    # return the results.
    task_hook.set_progress(2 / 3.0, "Formating results.")
    task_hook.set_results({
        "network": {"nodes": [], "edges": []},
        "node_attributes":
            {
                "node_types": [],
                "is_seed": [],
                "scores": []
        },
        'gene_interaction_dataset': "dummy",
        'drug_interaction_dataset': "dummy",
    })
