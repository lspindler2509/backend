from tasks.util.custom_network import add_edges
from tasks.task_hook import TaskHook
import graph_tool as gt
import gseapy as gp
from drugstone.models import *
from drugstone.serializers import *
import os
from drugstone.util.query_db import (
    query_proteins_by_identifier,
)


def parse_pathway(geneset, pathway, filtered_df, parameters, data_directory,background_mapping, background_mapping_reverse, map_genesets, gene_sets_dict, g = None):
    if isinstance(parameters, dict):
        id_space = parameters["config"].get("identifier", "symbol")
    else:
        parameters = json.loads(parameters)
        id_space = parameters["config"].get("identifier", "symbol")

    identifier_key = id_space
    if id_space == "ncbi":
        identifier_key = "entrez"
    elif id_space == "ensembl":
        identifier_key = "ensg"
    if g is None:
        ppi_dataset = parameters.get("ppi_dataset")
        pdi_dataset = parameters.get("pdi_dataset")
        custom_edges = parameters.get("custom_edges", False)
        filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
        if ppi_dataset['licenced'] or pdi_dataset['licenced']:
            filename += "_licenced"
        filename = os.path.join(data_directory, filename + ".gt")
        g = gt.load_graph(filename)
        if custom_edges:
            edges = parameters.get("input_network")['edges']
            g = add_edges(g, edges)
        

    
    seeds = list(parameters["seeds"])

    row = filtered_df.loc[(filtered_df['Gene_set'] == geneset) & (filtered_df['Term'] == pathway)].iloc[0]
    pathway = row['Term']
    genes = row['Genes']
    geneset = map_genesets[row['Gene_set']]
    genes = genes.split(";")
    only_pathway = list(set(gene_sets_dict[geneset][pathway]) - set(genes))
    filtered_only_pathway = []
    for gene in only_pathway:
        if gene in background_mapping:
            filtered_only_pathway.append(gene)
    only_pathway = filtered_only_pathway
    only_network = list(set(seeds) - set(genes))
    all_nodes = list(set(genes + only_pathway + only_network))
    nodes_mapped, identifier = query_proteins_by_identifier(all_nodes, identifier_key)
    nodes_mapped_dict = {node[identifier][0]: node for node in nodes_mapped}
    
    all_nodes_mapped = []
    isSeed = {}
    for node in all_nodes:
        drugstone_id = nodes_mapped_dict[node]["drugstone_id"]
        uniprot = nodes_mapped_dict[node]["uniprot"]
        symbol = nodes_mapped_dict[node]["symbol"]
        protein_name = nodes_mapped_dict[node]["protein_name"]
        entrez = nodes_mapped_dict[node]["entrez"]
        cellular_component = nodes_mapped_dict[node].get("cellular_component", [])
        layer = nodes_mapped_dict[node].get("layer", "")
        ensg = nodes_mapped_dict[node].get("ensg", "")
        if node in set(genes):
            isSeed[node] = True
            group = "overlap"
        elif node in set(only_pathway):
            isSeed[node] = False
            group = "onlyPathway"
        elif node in set(only_network):
            isSeed[node] = True
            group = "onlyNetwork"
            
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
    all_nodes_int = [int(background_mapping[gene]) for gene in all_nodes if gene in background_mapping]
    edges_unique = set()
    for node in all_nodes:
        for neighbor in g.get_all_neighbors(background_mapping[node]):
            if int(neighbor) > int(background_mapping[node]) and int(neighbor) in all_nodes_int:
                first_key = next(iter(background_mapping_reverse))

                if isinstance(first_key, int):
                    neighbor_key = int(neighbor)
                else:
                    neighbor_key = str(int(neighbor))
                edges_unique.add((node, background_mapping_reverse[neighbor_key]))
             
    final_network = {"nodes": all_nodes_mapped, "edges": [{"from": source, "to":target} for
                          source, target in edges_unique]}
    return final_network, isSeed



def add_group_to_config(config):
    if not config["node_groups"].get("overlap"):
        config["node_groups"]["overlap"] = {
                    "group_name": "overlap",
                    "color": {
                        "border": "#FC6C85",
                        "background": "#FC6C85",
                        "highlight": {
                            "border": "#FC6C85",
                            "background": "#FC6C85"
                        }
                    },
                    "shape": "circle",
                    "type": "gene",
                    "border_width": 0,
                    "border_width_selected": 0,
                    "font": {
                        "color": "#000000",
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
                    "group_id": "overlap"
                }
    if not config["node_groups"].get("onlyNetwork"):
        config["node_groups"]["onlyNetwork"] = {
                    "group_name": "only in network",
                    "color": {
                        "border": "#FFFF00",
                        "background": "#FFFF00",
                        "highlight": {
                            "border": "#FFFF00",
                            "background": "#FFFF00"
                        }
                    },
                    "shape": "circle",
                    "type": "gene",
                    "font": {
                        "color": "#000000",
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
                    "border_width": 1,
                    "border_width_selected": 2,
                    "shadow": True,
                    "group_id": "only_network"
                }
    if not config["node_groups"].get("onlyPathway"):
        config["node_groups"]["onlyPathway"]  = {
                    "group_name": "only in pathway",
                    "color": {
                        "border": "#FFCC09",
                        "background": "#FFCC09",
                        "highlight": {
                            "border": "#FFCC09",
                            "background": "#FFCC09"
                        }
                    },
                    "shape": "circle",
                    "type": "gene",
                    "font": {
                        "color": "#000000",
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
                    "border_width": 1,
                    "border_width_selected": 2,
                    "shadow": True,
                    "group_id": "only_pathway"
                }
        if not config["node_groups"].get("addedNode"):
            config["node_groups"]["addedNode"]  = {
                    "group_name": "added node",
                    "color": {
                        "border": "#FFB6C1",
                        "background": "#FFB6C1",
                        "highlight": {
                            "border": "#FFB6C1",
                            "background": "#FFB6C1"
                        }
                    },
                    "shape": "circle",
                    "type": "gene",
                    "font": {
                        "color": "#000000",
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
                    "border_width": 1,
                    "border_width_selected": 2,
                    "shadow": True,
                    "group_id": "added_node"
                }
    return config

def pathway_enrichment(task_hook: TaskHook):
    r"""
    Perform pathway enrichment analysis.

    Parameters
    ----------

    num_threads : int, optional (default: 1)
      Number of threads. Requires that graph_tool is compiled with OpenMP support.
      Should not be exposed in the frontend.
      
    alpha : number, optional (default: 0.05)
    The cutoff for the adjusted p-value of the enriched pathways.
    
    kegg: bool, optional (default: True)
    If True, KEGG pathways are included in the analysis.
    
    reactome: bool, optional (default: True)
    If True, Reactome pathways are included in the analysis.
    
    wiki: bool, optional (default: True)
    If True, Wiki pathways are included in the analysis.
    
      
    Returns
    -------
    results : {
        "algorithm": "pathway_enrichment",
        "filteredDf": filtered_df.to_json(orient='records'),
        "backgroundMapping": background_mapping,
        "backgroundMappingReverse": background_mapping_reverse,
        "mapGenesets": map_genesets,
        "mapGenesetsReverse": map_genesets_reverse,
        "geneSetsDict": gene_sets_dict,
        "table_view": table_view_results,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
        "parameters": task_hook.parameters,
        "geneSetPathways": gene_set_terms_dict,
        "config": add_group_to_config(task_hook.parameters["config"]),
    }
    
    "algorithm": "pathway_enrichment"
    "filteredDf": A JSON string containing the filtered DataFrame with the pathway enrichment results.
    "backgroundMapping": A dictionary that maps the internal node IDs to the original node IDs.
    "backgroundMappingReverse": A dictionary that maps the original node IDs to the internal node IDs.
    "mapGenesets": A dictionary that maps the internal geneset IDs to the original geneset IDs.
    "mapGenesetsReverse": A dictionary that maps the original geneset IDs to the internal geneset IDs.
    "geneSetsDict": A dictionary that contains the genesets and their pathways.
    "table_view": A list of dictionaries containing the enriched pathways.
    "gene_interaction_dataset": The gene interaction dataset.
    "drug_interaction_dataset": The drug interaction dataset.
    "parameters": The parameters of the task.
    "geneSetPathways": A dictionary that contains the genesets and their pathways.
    "config": The configuration of the task.

    Notes
    -----
    This implementation is based on graph-tool, a very efficient Python package for network
    analysis with C++ backend and multi-threading support. Installation instructions for graph-tool
    can be found at https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions.

    References
    ----------
    ..  [1] Fang, Zhuoqing, Xinyuan Liu, and Gary Peltz. “GSEApy: A Comprehensive Package for Performing Gene Set Enrichment Analysis in Python.” Bioinformatics (Oxford, England) 39, no. 1 (January 1, 2023): btac757. https://doi.org/10.1093/bioinformatics/btac757.
        [2] Xie, Zhuorui, Allison Bailey, Maxim V. Kuleshov, Daniel J. B. Clarke, John E. Evangelista, Sherry L. Jenkins, Alexander Lachmann, et al. “Gene Set Knowledge Discovery with Enrichr.” Current Protocols 1, no. 3 (March 2021): e90. https://doi.org/10.1002/cpz1.90.
    """
    
    # Get the data directory.
    data_dir = os.path.dirname(os.path.dirname(task_hook.data_directory))

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
    
    # Type: number.
    # Semantics: Alpha value as cutoff for the adjusted p-value.
    # Example: 0.05
    # Reasonable default: 0.05
    # Note: The user can specify this value when starting the task.
    alpha = task_hook.parameters.get("alpha", 0.05)

    # Parsing input file.
    task_hook.set_progress(1 / 4.0, "Parsing input.")
    
    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    identifier_key = id_space
    if id_space == "ncbi":
        identifier_key = "entrez"
    elif id_space == "ensembl":
        identifier_key = "ensg"
    
    filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset['licenced'] or pdi_dataset['licenced']:
        filename += "_licenced"
    filename = os.path.join(task_hook.data_directory, filename + ".gt")
    g = gt.load_graph(filename)
    if custom_edges:
        edges = task_hook.parameters.get("input_network")['edges']
        g = add_edges(g, edges)
    
    background = []
    background_mapping = {}
    background_mapping_reverse = {}
    for v in g.vertices():
        if g.vertex_properties["type"][v] == "protein":
            background.append(g.vertex_properties["internal_id"][v])
            background_mapping[g.vertex_properties["internal_id"][int(v)]] = int(v)
            background_mapping_reverse[int(v)] = g.vertex_properties["internal_id"][(v)]
    
    gene_sets = []
    gene_sets_dict = {}
    map_genesets = {}
    map_genesets_reverse = {}
    
    # parse genesets
    if task_hook.parameters.get("kegg", True):
        pathway_kegg = {}
        path = os.path.join(data_dir, "gene_sets", "kegg_"+identifier_key+".txt")
        with open(path, "r") as file:
            for line in file:
                parts = line.strip().split('\t')
                pathway = parts[0]
                genes = parts[1:]
                genes = [gene for gene in genes if gene]
                pathway_kegg[pathway] = genes
        gene_sets.append(pathway_kegg)
        gene_sets_dict["kegg"] = pathway_kegg
        map_genesets["gs_ind_0"] = "kegg"
        map_genesets_reverse["kegg"] = "gs_ind_0"
    
    if task_hook.parameters.get("reactome", True):
        pathway_reactome = {}
        path = os.path.join(data_dir, "gene_sets", "reactome_"+identifier_key+".txt")
        with open(path, "r") as file:
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
            map_genesets_reverse["reactome"] = "gs_ind_1"
        else:
            map_genesets["gs_ind_0"] = "reactome"
            map_genesets_reverse["reactome"] = "gs_ind_0"
        

    
    if task_hook.parameters.get("wiki", True):
        pathway_wiki = {}
        path = os.path.join(data_dir, "gene_sets", "wiki_"+identifier_key+".txt")
        with open(path, "r") as file:
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
            map_genesets_reverse["wiki"] = "gs_ind_0"
        elif len(map_genesets.keys()) == 1:
            map_genesets["gs_ind_1"] = "wiki"
            map_genesets_reverse["wiki"] = "gs_ind_1"
        else:
            map_genesets["gs_ind_2"] = "wiki"
            map_genesets_reverse["wiki"] = "gs_ind_2"

    task_hook.set_progress(2 / 4.0, "Running pathway enrichment.")

    # Perform pathway enrichment analysis.
    enr = gp.enrichr(gene_list=seeds,
                     gene_sets=gene_sets,
                     organism='human',
                     outdir=None,
                     background=background,
                     )
    
    task_hook.set_progress(3 / 4.0, "Parse pathway enrichment result for lowest adjusted p-value.")
          
    # filter result accroding to adjusted p-value and sort by adjusted p-value
    filtered_df = enr.results[enr.results['Adjusted P-value'] <= alpha]
    filtered_df = filtered_df.sort_values(by=['Adjusted P-value'])
          
    # parse data for tableview      
    table_view_results = []
    for _ , row in filtered_df.iterrows():
        geneset = map_genesets[row['Gene_set']]
        pathway = row['Term']
        table_view_results.append({"geneset": geneset, "pathway": pathway, "overlap": row['Overlap'], "adj_pvalue": row['Adjusted P-value'], "odds_ratio": round(row['Odds Ratio'], 2), "genes": row['Genes']})

    gene_sets_list = filtered_df['Gene_set'].unique().tolist()
    gene_set_terms_dict = {}
    genesets = []
    for gene_set in gene_sets_list:
        geneset = map_genesets[gene_set]
        genesets.append(geneset)
        terms_list = filtered_df[filtered_df['Gene_set'] == gene_set]['Term'].tolist()
        gene_set_terms_dict[geneset] = terms_list
    
   
    # return the results.
    task_hook.set_progress(4 / 4.0, "Formating results.")
        
    # keep calculations that can be used again when user chooses different pathway
    task_hook.set_results({
        "algorithm": "pathway_enrichment",
        "filteredDf": filtered_df.to_json(orient='records'),
        "backgroundMapping": background_mapping,
        "backgroundMappingReverse": background_mapping_reverse,
        "mapGenesets": map_genesets,
        "mapGenesetsReverse": map_genesets_reverse,
        "geneSetsDict": gene_sets_dict,
        "table_view": table_view_results,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
        "parameters": task_hook.parameters,
        "geneSetPathways": gene_set_terms_dict,
        "config": add_group_to_config(task_hook.parameters["config"]),
    })
