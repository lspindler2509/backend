import graph_tool as gt
import graph_tool.topology as gtt


# def read_graph_tool_graph(file_path, seeds, datasets, ignored_edge_types, max_deg, ignore_non_seed_baits=False, include_indirect_drugs=False, include_non_approved_drugs=False):
def read_graph_tool_graph(file_path, seeds, id_space, max_deg, include_indirect_drugs=False,
                          include_non_approved_drugs=False,
                          target='drug'):
    r"""Reads a graph-tool graph from file.

    Reads a graph-tool graph from graphml or gt file and returns is along
    with the internal IDs of the seed and viral seeds and the drugs.

    Parameters
    ----------
    file_path : str | graph-tool Graph
      A string specifying the path to a graphml or gt file, or graph tool graph.

    seeds : list of str
      A list of drugstone IDs identifying the seed seeds.

    include_indirect_drugs : bool
      If True, edges from non-seed host proteins to drugs are ignored when ranking drugs.

    include_non_approved_drugs : bool
      If True, also non-approved drugs are included in the analysis

    target : str
      A string specifying the target of the search, either "drug" or "drug-target"

    Returns
    -------
    g : graph_tool.Graph
      The constructed graph.

    seed_ids : list of int
      The graph indices for all seed nodes

    drug_ids : list of int
      The graph indices for all drug nodes
    """
    # Read the graph.
    if isinstance(file_path, str):
        g = gt.load_graph(file_path)
    else:
        g = file_path
        
    # drug_protein = "DrugHasTarget"
    d_type = "drug"
    node_name_attribute = "internal_id"  # nodes in the input network which is created from RepoTrialDB have primaryDomainId as name attribute
    # Delete all nodes that are not contained in the selected datasets and have degrees higher than max_deg
    deleted_nodes = []
    for node in range(g.num_vertices()):
        # Remove all unconnected nodes TODO probably already skip when creating .gt files
        if g.vertex(node).out_degree() == 0 and target == 'drug':
            deleted_nodes.append(node)
        elif not g.vertex_properties[node_name_attribute][node] in set(seeds) and (
                g.vertex(node).out_degree() > max_deg):
            deleted_nodes.append(node)
        # remove all drugs from graph if we are not looking for drugs
        elif target != 'drug' and g.vertex_properties["type"][node] == d_type:
            deleted_nodes.append(node)

    g.remove_vertex(reversed(sorted(deleted_nodes)), fast=True)

    # Retrieve internal IDs of seed_ids
    seeds = set(seeds)
    seed_ids = {}
    drug_ids = []
    for node in range(g.num_vertices()):
        node_type = g.vertex_properties["type"][node]
        seed_id = g.vertex_properties[node_name_attribute][node]
        if seed_id in seeds:
            seed_ids[node] = seed_id
        if node_type == d_type:
            if include_non_approved_drugs:
                drug_ids.append(node)
            else:
                # drug_groups = g.vertex_properties["status"][node].split(', ')
                if "approved" in g.vertex_properties["status"][node]:
                    drug_ids.append(node)

    # Delete edges that should be ignored or are not contained in the selected dataset.
    deleted_edges = []

    for edge in g.edges():
        if edge.source == edge.target:
            deleted_edges.append(edge)

    g.set_fast_edge_removal(fast=True)
    for edge in deleted_edges:
        g.remove_edge(edge)
    g.set_fast_edge_removal(fast=False)

    deleted_edges = []

    # If only_direct_drugs should be included, remove any drug-protein edges that the drug is not a direct neighbor of
    # any seeds
    if drug_ids and not include_indirect_drugs:
        direct_drugs = set()
        for edge in g.edges():
            if g.vertex_properties["type"][edge.target()] == d_type and edge.source() in seed_ids:
                direct_drugs.add(edge.target())
            elif g.vertex_properties["type"][edge.source()] == d_type and edge.target() in seed_ids:
                direct_drugs.add(edge.source())
        for edge in g.edges():
            if g.edge_properties["type"][edge] == 'drug-protein':
                if g.vertex_properties["type"][edge.target()] == d_type:
                    indir_drug = edge.target() not in direct_drugs
                    not_seed = edge.source() not in seed_ids
                    if indir_drug or not_seed:
                        deleted_edges.append(edge)
                    if indir_drug and int(edge.target()) in drug_ids:
                        drug_ids.remove(int(edge.target()))

                elif g.vertex_properties["type"][edge.source()] == d_type and \
                        edge.source() not in direct_drugs or edge.target() not in seed_ids:
                    indir_drug = edge.source() not in direct_drugs
                    not_seed = edge.target() not in seed_ids
                    if indir_drug or not_seed:
                        deleted_edges.append(edge)
                    if indir_drug and int(edge.source()) in drug_ids:
                        drug_ids.remove(int(edge.source()))
            # else:
            #     deleted_edges.append(edge)

    g.set_fast_edge_removal(fast=True)
    for edge in deleted_edges:
        g.remove_edge(edge)
    g.set_fast_edge_removal(fast=False)

    # Return the graph and the indices of the seed_ids and the seeds.
    return g, list(seed_ids.keys()), drug_ids
