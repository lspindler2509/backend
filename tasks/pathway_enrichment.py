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
import os


def pathway_enrichment(task_hook: TaskHook):
    r"""
    Perform pathway enrichment analysis.

    Parameters
    ----------
    

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
    data_dir = os.path.dirname(os.path.dirname(task_hook.data_directory))
    PATH_REACTOME_GENESET = os.path.join(data_dir, "gene_sets/KEGG_2021_Human.txt")
    PATH_WIKI_GENESET = os.path.join(data_dir, "gene_sets/Reactome_2022.txt")
    PATH_KEGG_GENESET = os.path.join(data_dir, "gene_sets/WikiPathway_2023_Human.txt")
    
    node_name_attribute = "internal_id" # nodes in the input network which is created from RepoTrialDB have primaryDomainId as name attribute


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

    search_target = task_hook.parameters.get("target", "gene")

    id_space = task_hook.parameters["config"].get("identifier", "symbol")

    custom_edges = task_hook.parameters.get("custom_edges", False)
    
    alpha = task_hook.parameters.get("alpha", 0.05)

    # Parsing input file.
    task_hook.set_progress(0 / 3.0, "Parsing input.")
    
    # we always want to use symbol as id space for the enrichment analysis
    filename = f"symbol_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset['licenced'] or pdi_dataset['licenced']:
        filename += "_licenced"
    filename = os.path.join(task_hook.data_directory, filename + ".gt")
    g = gt.load_graph(filename)
    
    background = []
    background_mapping = {}
    background_mapping_reverse = {}
    for v in g.vertices():
        if g.vertex_properties["type"][v] == "protein":
            background.append(g.vertex_properties["internal_id"][v])
            background_mapping[g.vertex_properties["internal_id"][int(v)]] = int(v)
            background_mapping_reverse[int(v)] = g.vertex_properties["internal_id"][(v)]


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
            if id_space == "entrez" or id_space == "ncbigene":
                seed_without_identifier = seed.replace("entrez.", "")
                protein = models.Protein.objects.filter(
                    entrez=seed_without_identifier
                ).last()
                seeds_symbol.append(protein.gene)
            if id_space == "ensembl" or id_space == "ensg":
                print("Not parsed yet! Map via EnsemblGene")

    print(seeds_symbol)
    
    gene_sets = []
    gene_sets_dict = {}
    map_genesets = {}
    
    # parse genesets
    
    if task_hook.parameters.get("kegg"):
        pathway_kegg = {}
        with open(PATH_KEGG_GENESET, "r") as file:
            for line in file:
                parts = line.strip().split('\t')
                pathway = parts[0]
                genes = parts[1:]
                genes = [gene for gene in genes if gene]
                pathway_kegg[pathway] = genes
        gene_sets.append(pathway_kegg)
        gene_sets_dict["kegg"] = pathway_kegg
        map_genesets["gs_ind_0"] = "kegg"
    
    if task_hook.parameters.get("reactome"):
        pathway_reactome = {}
        with open(PATH_REACTOME_GENESET, "r") as file:
            for line in file:
                parts = line.strip().split('\t')
                pathway = parts[0]
                genes = parts[1:]
                genes = [gene for gene in genes if gene]
                pathway_reactome[pathway] = genes
        gene_sets.append(pathway_reactome)
        gene_sets_dict["reactome"] = pathway_reactome
        if map_genesets.get("gs_ind_0", False):
            map_genesets["gs_ind_1"] = "reactome"
        else:
            map_genesets["gs_ind_0"] = "reactome"
        

    
    if task_hook.parameters.get("wiki"):
        pathway_wiki = {}
        with open(PATH_WIKI_GENESET, "r") as file:
            for line in file:
                parts = line.strip().split('\t')
                pathway = parts[0]
                genes = parts[1:]
                genes = [gene for gene in genes if gene]
                pathway_wiki[pathway] = genes
        gene_sets.append(pathway_wiki)
        gene_sets_dict["wiki"] = pathway_wiki
        
        if len(map_genesets.keys()) == 0:
            map_genesets["gs_ind_0"] = "wiki"
        elif len(map_genesets.keys()) == 1:
            map_genesets["gs_ind_1"] = "wiki"
        else:
            map_genesets["gs_ind_2"] = "wiki"

        

    enr = gp.enrichr(gene_list=seeds_symbol,
                     gene_sets=gene_sets,
                     organism='human',
                     outdir=None,
                     background=background,
                     )
    
    # filter result accroding to adjusted p-value
    filtered_df = enr.results[enr.results['Adjusted P-value'] <= alpha]
    
    networks = {}
    
    # 3 groups: overlap, only_network, only_pathway
    
    former_network = task_hook.parameters.get("input_network")
    
    genes_not_in_network = set()
    
    for index, row in filtered_df.iterrows():
        pathway = row['Term']
        genes = row['Genes']
        geneset = map_genesets[row['Gene_set']]
        genes = genes.split(";")
        only_pathway = list(set(gene_sets_dict[geneset][pathway]) - set(genes))
        filtered_only_pathway = []
        for gene in only_pathway:
            if gene in background_mapping:
                genes_not_in_network.add(gene)
                # TODO: is that what we want?
                filtered_only_pathway.append(gene)
        only_pathway = filtered_only_pathway
        only_network = []
        for node in former_network["nodes"]:
            only_network.extend(node.get("symbol", []))
        only_network = list(set(only_network) - set(genes))
        all_nodes = list(set(genes + only_pathway + only_network))
        all_nodes_int = [int(background_mapping[gene]) for gene in all_nodes if gene in background_mapping]
        edges_unique = set()
        for node in all_nodes:
            for neighbor in g.get_all_neighbors(background_mapping[node]):
                if int(neighbor) > int(background_mapping[node]) and int(neighbor) in all_nodes_int:
                    edges_unique.add((node, background_mapping_reverse[int(neighbor)]))
             
        
        networks.setdefault(geneset, {})[pathway] = {"nodes": all_nodes, "edges": [{"from": source, "to":target} for
                          source, target in edges_unique]}

    print(networks)
    
    # extract edges from the network
    
    
    # TODO: namespace????
    if custom_edges:
        edges = task_hook.parameters.get("input_network")['edges']
        g = add_edges(g, edges)
        
        
        
    print(len(filtered_df))

    # return the results.
    task_hook.set_progress(2 / 3.0, "Formating results.")
    task_hook.set_results({
        "network": {"nodes": [], "edges": []},
        "node_attributes":
            {
                "node_types": [],
                "is_seed": [],
                "scores": [],
                "details": []
        },
        'gene_interaction_dataset': "dummy",
        'drug_interaction_dataset': "dummy",
    })
