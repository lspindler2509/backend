from collections import defaultdict
from typing import List, Tuple
import graph_tool.all as gt
from drugstone import models
import multiprocessing
from django import db
from pathlib import Path
from django.core.management import BaseCommand
import django
import os

django.setup()

KERNEL = int(os.environ.get('GT_THREADS', 6))


def _internal_expression_scores(drugstone_id: str) -> dict:
    """ Looks up the tissue specific expression scores for a given protein.
    The scores are loaded from the django database.

    Args:
        drugstone_id (str): drugstone id of protein in format 'pxxxx'

    Returns:
        dict: keys are tissue-names and values are the expression scores
    """
    protein_object = models.Protein.objects.get(id=int(drugstone_id[1:]))

    # get expression scores
    tissues = models.Tissue.objects.all()
    tissue_scores = {t.name: None for t in tissues}
    for t in tissues:
        res = models.ExpressionLevel.objects.filter(
            tissue=t,
            protein=protein_object
        )
        if res:
            tissue_scores[t.name] = res[0].expression_level

    return tissue_scores


def _internal_pdi(dataset) -> List[models.ProteinDrugInteraction]:
    """ Fetches all internal protein-drug interactions for a given dataset.
    Interactions are taken from the django database.

    Args:
        dataset_name (str): Name of the dataset, e.g. "DrugBank"

    Returns:
        List[dict]: List of representaions of interaction objects
    """
    # get all interactions
    node_node_interaction_objects = models.ProteinDrugInteraction.objects.filter(
        pdi_dataset__id=dataset.id
    )
    # node_node_interactions = serializers.ProteinDrugInteractionSerializer(many=True) \
    #     .to_representation(node_node_interaction_objects)

    return node_node_interaction_objects

def _internal_pdis(dataset) -> List[models.ProteinDisorderAssociation]:
    """ Fetches all internal protein-disorder associations for a given dataset.
    Interactions are taken from the django database.

    Args:
        dataset_name (str): Name of the dataset, e.g. "DrugBank"

    Returns:
        List[dict]: List of representaions of interaction objects
    """
    # get all interactions
    node_node_interaction_objects = models.ProteinDisorderAssociation.objects.filter(
        pdis_dataset__id=dataset.id
    )

    return node_node_interaction_objects


def _internal_drdis(dataset) -> List[models.DrugDisorderIndication]:
    node_node_interaction_objects = models.DrugDisorderIndication.objects.filter(
        drdi_dataset__id=dataset.id
    )

    return node_node_interaction_objects


def _internal_ppis(dataset) -> List[models.ProteinProteinInteraction]:
    """ Fetches all internal protein-protein interactions for a given dataset.
    Interactions are taken from the django database.

    Args:
        dataset_name (str): Name of the dataset, e.g. "BioGRID"

    Returns:
        List[dict]: List of representaions of interaction objects
    """
    # get all interactions
    node_node_interaction_objects = models.ProteinProteinInteraction.objects.filter(
        ppi_dataset__id=dataset.id
    )

    return node_node_interaction_objects

def get_filename(dataset_name, dataset_version, identifier=None, edge_type="", licensed=False, isReviewed=False, fmt="gt"):
    filename = f"./data/Networks/{dataset_name}_{dataset_version}_"

    if isReviewed:
        filename += "reviewed-"

    if identifier is not None:
        filename += identifier+"-"
    filename += edge_type

    if licensed:
        filename += "_licenced"

    filename += "_download." + fmt
    return filename


def get_or_create_pdi_network(dataset, identifier, licensed, isReviewed, fmt):

    # dataset, dataset_type, identifier, isReviewed, licensed, format = params


    filename = get_filename(dataset_name=dataset.name, dataset_version=dataset.version, identifier=identifier, edge_type="protein-drug-interaction", licensed=licensed, fmt=fmt, isReviewed=isReviewed)

    filepath = Path(filename)
    if os.path.exists(filepath):
        return filepath

    print(f'Creating {filename}')

    g = gt.Graph(directed=False)

    # For edges
    e_type = g.new_edge_property("string")
    g.edge_properties["type"] = e_type

    e_actions=g.new_edge_property("string")
    g.edge_properties["actions"] = e_actions

    # for all nodes

    v_type = g.new_vertex_property("string")
    g.vertex_properties["type"] = v_type

    v_name = g.new_vertex_property("string")
    g.vertex_properties["label"] = v_name

    v_internal_id = g.new_vertex_property("string")
    g.vertex_properties["internal_ids"] = v_internal_id

    v_id = g.new_vertex_property("string")
    g.vertex_properties["id"] = v_id

    # For protein nodes
    v_reviewed = g.new_vertex_property("boolean")
    g.vertex_properties["reviewed"] = v_reviewed

    # For drug nodes
    v_status = g.new_vertex_property("string")
    g.vertex_properties["status"]  = v_status
    # store nodes to connect them when creating edges
    vertices = {}
    drug_vertices = {}
    # add vertices

    print(f'loading nodes for {identifier}')

    is_entrez = (identifier == 'entrez' or identifier == 'ncbigene')
    is_symbol = identifier == 'symbol'
    is_uniprot = identifier == 'uniprot'
    is_ensg = (identifier == 'ensg' or identifier == 'ensembl')
    is_internal = identifier == None

    if is_ensg or is_internal:
        ensembl_set = defaultdict(set)
        for node in models.EnsemblGene.objects.all():
            ensembl_set[node.protein_id].add(node.name)

    node_id_map = defaultdict(set)
    drugstone_id_to_node = dict()

    if isReviewed:
        proteins = models.Protein.objects.filter(isReviewed=True)
    else:
        proteins = models.Protein.objects.all()

    for node in proteins:
        if is_entrez:
            if len(node.entrez) != 0:
                node_id_map[node.entrez].add(node.id)
                drugstone_id_to_node[node.id] = node
        elif is_symbol:
            if len(node.gene) != 0:
                node_id_map[node.gene].add(node.id)
                drugstone_id_to_node[node.id] = node
        elif is_uniprot:
            node_id_map[node.uniprot_code].add(node.id)
            drugstone_id_to_node[node.id] = node
        elif is_ensg:
            for id in ensembl_set[node.id]:
                node_id_map[id].add(node.id)
                drugstone_id_to_node[node.id] = node
        elif is_internal:
            node_id_map[node.id].add(node.id)
            drugstone_id_to_node[node.id] = node

            v_uniprot = g.new_vertex_property("string")
            g.vertex_properties["uniprot"] = v_uniprot

            v_symbol = g.new_vertex_property("string")
            g.vertex_properties["symbol"] = v_symbol

            v_entrez = g.new_vertex_property("string")
            g.vertex_properties["entrez"] = v_entrez

            v_ensembl = g.new_vertex_property("string")
            g.vertex_properties["ensembl"] = v_ensembl
        print(f"Protein nodes:{len(node_id_map.items())}")
        done = 0
        for id, internal_ids in node_id_map.items():
            print(f"protein: {id}")
            v = g.add_vertex()
            v_type[v] = 'protein'
            v_internal_id[v] = ",".join({f"p{id}" for id in internal_ids})
            if is_internal:
                node = drugstone_id_to_node[id]
                v_reviewed[v] = node.isReviewed
                v_name[v] = node.gene
                v_uniprot[v] = node.uniprot_code
                v_symbol[v] = node.gene
                v_entrez[v] = node.entrez
                if node.id in ensembl_set.keys():
                    v_ensembl[v] = ",".join({f"{id}" for id in ensembl_set[node.id]})
                done+=1
                print(f"done: {done}/{len(node_id_map.items())}")
                vertices[id] = v
            else:
                for drugstone_id in internal_ids:
                    print(f"\tfor {drugstone_id}")
                    node = drugstone_id_to_node[drugstone_id]
                    v_reviewed[v] = node.isReviewed
                    v_name[v] = node.gene
                    vertices[drugstone_id] = v


    for node in models.Drug.objects.all():
        v = g.add_vertex()
        v_type[v] = 'drug'
        v_name[v] = node.name
        v_status[v] = node.status
        v_internal_id[v] = f'dr{node.id}'

        drug_vertices[node.id] = v

    uniq_edges = set()
    n = 0
    print(f'loading drug_edges/{dataset}')
    for edge_raw in _internal_pdi(dataset):
        id1 = edge_raw.drug_id
        id2 = edge_raw.protein_id
        hash = f'{id1}_{id2}'
        if hash not in uniq_edges and id1 in drug_vertices and id2 in vertices:
            uniq_edges.add(hash)
            e = g.add_edge(drug_vertices[id1], vertices[id2])
            n += 1
            e_type[e] = 'drug-protein'
            e_actions[e] = edge_raw.actions
    print("done with drug edges: ", n)

    # remove unconnected proteins
    delete_vertices = set()
    for vertex in vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)
    print("removing unconnected proteins: ", len(delete_vertices))

    # remove unconnected drugs
    for vertex in drug_vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)

    g.remove_vertex(reversed(sorted(delete_vertices)), fast=True)
    Path('./data/Networks/').mkdir(parents=True, exist_ok=True)
    g.save(filename, fmt=fmt)
    print(f"Created file {filepath}")
    print("Size of graph - nodes: ", g.num_vertices(), " edges: ", g.num_edges())
    return filepath

def get_or_create_drdis_network(dataset, licensed, fmt):


    filename = get_filename(dataset_name=dataset.name, dataset_version=dataset.version, edge_type="drug-disorder-indication", licensed=licensed, fmt=fmt)


    filepath = Path(filename)
    if os.path.exists(filepath):
        return filepath

    print(f'Creating {filename}')

    g = gt.Graph(directed=False)

    # For edges
    e_type = g.new_edge_property("string")
    g.edge_properties["type"] = e_type


    # for all nodes
    v_type = g.new_vertex_property("string")
    g.vertex_properties["type"] = v_type

    v_name = g.new_vertex_property("string")
    g.vertex_properties["label"] = v_name

    v_internal_id = g.new_vertex_property("string")
    g.vertex_properties["internal_ids"] = v_internal_id

    v_id = g.new_vertex_property("string")
    g.vertex_properties["id"] = v_id

    # for disorder nodes
    v_icd10 = g.new_vertex_property("string")
    g.vertex_properties["icd10_code"] = v_icd10

    # For drug nodes
    v_status = g.new_vertex_property("string")
    g.vertex_properties["status"]  = v_status
    # store nodes to connect them when creating edges
    disorder_vertices = {}
    drug_vertices = {}
    # add vertices


    for node in models.Disorder.objects.all():
        v = g.add_vertex()
        v_type[v] = 'disorder'
        v_name[v] = node.label
        v_icd10[v] = node.icd10
        v_internal_id[v] = f'dis{node.id}'

        disorder_vertices[node.id] = v

    for node in models.Drug.objects.all():
        v = g.add_vertex()
        v_type[v] = 'drug'
        v_name[v] = node.name
        v_status[v] = node.status
        v_internal_id[v] = f'dr{node.id}'

        drug_vertices[node.id] = v

    uniq_edges = set()
    n = 0
    print(f'loading drug-disease_edges/{dataset}')
    for edge_raw in _internal_drdis(dataset):
        id1 = edge_raw.drug_id
        id2 = edge_raw.disorder_id
        hash = f'{id1}_{id2}'
        if hash not in uniq_edges and id1 in drug_vertices and id2 in disorder_vertices:
            uniq_edges.add(hash)
            e = g.add_edge(drug_vertices[id1], disorder_vertices[id2])
            n += 1
            e_type[e] = 'drug-disorder'
    print("done with drug-disorder edges: ", n)

    # remove unconnected proteins
    delete_vertices = set()
    for vertex in disorder_vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)
    print("removing unconnected disorder: ", len(delete_vertices))

    # remove unconnected drugs
    for vertex in drug_vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)
    print("removing unconnected drugs: ", len(delete_vertices))

    g.remove_vertex(reversed(sorted(delete_vertices)), fast=True)
    Path('./data/Networks/').mkdir(parents=True, exist_ok=True)
    g.save(filename, fmt=fmt)
    print(f"Created file {filepath}")
    print("Size of graph - nodes: ", g.num_vertices(), " edges: ", g.num_edges())
    return filepath

def get_or_create_pdis_network(dataset, identifier, licensed, isReviewed, fmt):

    # dataset, dataset_type, identifier, isReviewed, licensed, format = params

    filename = get_filename(dataset_name=dataset.name, dataset_version=dataset.version, identifier=identifier,
                            edge_type="protein-disorder-association", licensed=licensed, fmt=fmt, isReviewed=isReviewed)

    filepath = Path(filename)
    if os.path.exists(filepath):
        return filepath

    print(f'Creating {filename}')

    g = gt.Graph(directed=False)

    #For edges
    e_type = g.new_edge_property("string")
    g.edge_properties["type"] = e_type

    e_score = g.new_edge_property("float")
    g.edge_properties["score"] = e_score


#for all nodes

    v_type = g.new_vertex_property("string")
    g.vertex_properties["type"] = v_type

    v_name = g.new_vertex_property("string")
    g.vertex_properties["label"] = v_name

    v_internal_id = g.new_vertex_property("string")
    g.vertex_properties["internal_ids"] = v_internal_id

    v_id = g.new_vertex_property("string")
    g.vertex_properties["id"] = v_id


    #For protein nodes
    v_reviewed = g.new_vertex_property("boolean")
    g.vertex_properties["reviewed"] = v_reviewed

    # for disorder nodes
    v_icd10 = g.new_vertex_property("string")
    g.vertex_properties["icd10_code"] = v_icd10


    # store nodes to connect them when creating edges
    vertices = {}
    disorder_vertices = {}
    # add vertices

    print(f'loading nodes for {identifier}')

    is_entrez = (identifier == 'entrez' or identifier == 'ncbigene')
    is_symbol = identifier == 'symbol'
    is_uniprot = identifier == 'uniprot'
    is_ensg = (identifier == 'ensg' or identifier == 'ensembl')
    is_internal = identifier == None

    if is_ensg or is_internal:
        ensembl_set = defaultdict(set)
        for node in models.EnsemblGene.objects.all():
            ensembl_set[node.protein_id].add(node.name)

    node_id_map = defaultdict(set)
    drugstone_id_to_node = dict()

    if isReviewed:
        proteins = models.Protein.objects.filter(isReviewed=True)
    else:
        proteins = models.Protein.objects.all()

    for node in proteins:
        if is_entrez:
            if len(node.entrez) != 0:
                node_id_map[node.entrez].add(node.id)
                drugstone_id_to_node[node.id] = node
        elif is_symbol:
            if len(node.gene) != 0:
                node_id_map[node.gene].add(node.id)
                drugstone_id_to_node[node.id] = node
        elif is_uniprot:
            node_id_map[node.uniprot_code].add(node.id)
            drugstone_id_to_node[node.id] = node
        elif is_ensg:
            for id in ensembl_set[node.id]:
                node_id_map[id].add(node.id)
                drugstone_id_to_node[node.id] = node
        elif is_internal:
            node_id_map[node.id].add(node.id)
            drugstone_id_to_node[node.id] = node

            v_uniprot = g.new_vertex_property("string")
            g.vertex_properties["uniprot"] = v_uniprot

            v_symbol = g.new_vertex_property("string")
            g.vertex_properties["symbol"] = v_symbol

            v_entrez = g.new_vertex_property("string")
            g.vertex_properties["entrez"] = v_entrez

            v_ensembl = g.new_vertex_property("string")
            g.vertex_properties["ensembl"] = v_ensembl



    for id, internal_ids in node_id_map.items():
        v = g.add_vertex()
        v_type[v] = 'protein'
        v_internal_id[v] = ",".join({f"p{id}" for id in internal_ids})
        for drugstone_id in internal_ids:
            node = drugstone_id_to_node[drugstone_id]
            v_reviewed[v] = node.isReviewed
            v_name[v] = node.gene
            vertices[drugstone_id] = v
        if is_internal:
            node = drugstone_id_to_node[id]
            v_uniprot[v] = node.uniprot_code
            v_symbol[v] = node.gene
            v_entrez[v] = node.entrez
            if node.id in ensembl_set.keys():
                v_ensembl[v]  = ",".join({f"{id}" for id in ensembl_set[node.id]})


    for node in models.Disorder.objects.all():
        v = g.add_vertex()
        v_type[v] = 'disorder'
        v_name[v] = node.label
        v_icd10[v] = node.icd10
        v_internal_id[v] = f'dis{node.id}'

        disorder_vertices[node.id] = v

    uniq_edges = set()
    n = 0
    print(f'loading protein-disorder_edges/{dataset}')
    for edge_raw in _internal_pdis(dataset):
        id1 = edge_raw.disorder_id
        id2 = edge_raw.protein_id
        hash = f'{id1}_{id2}'
        if hash not in uniq_edges and id1 in disorder_vertices and id2 in vertices:
            uniq_edges.add(hash)
            e = g.add_edge(disorder_vertices[id1], vertices[id2])
            n += 1
            e_type[e] = 'protein-disorder'
            e_score[e] = edge_raw.score
    print("done with protein-disorder edges: ", n)

    # remove unconnected proteins
    delete_vertices = set()
    for vertex in vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)
    print("removing unconnected proteins: ", len(delete_vertices))

    # remove unconnected drugs
    for vertex in disorder_vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)

    g.remove_vertex(reversed(sorted(delete_vertices)), fast=True)
    Path('./data/Networks/').mkdir(parents=True, exist_ok=True)
    g.save(filename, fmt=fmt)
    print(f"Created file {filepath}")
    print("Size of graph - nodes: ", g.num_vertices(), " edges: ", g.num_edges())
    return filepath

def get_or_create_ppi_network(dataset, identifier, licensed, isReviewed, fmt):

    filename = get_filename(dataset_name=dataset.name, dataset_version=dataset.version, identifier=identifier,
                            edge_type="protein-protein-interaction", licensed=licensed, fmt=fmt, isReviewed=isReviewed)

    filepath = Path(filename)
    if os.path.exists(filepath):
        return filepath

    print(f'Creating {filename}')

    g = gt.Graph(directed=False)

    # For edges
    e_type = g.new_edge_property("string")
    g.edge_properties["type"] = e_type

    e_directed = g.new_edge_property("bool")
    g.edge_properties["directed"] = e_directed

    e_stimulation = g.new_edge_property("bool")
    g.edge_properties["stimulation"] = e_stimulation

    e_inhibition = g.new_edge_property("bool")
    g.edge_properties["inhibition"] = e_inhibition


    #for all nodes

    v_type = g.new_vertex_property("string")
    g.vertex_properties["type"] = v_type

    v_name = g.new_vertex_property("string")
    g.vertex_properties["label"] = v_name

    v_internal_id = g.new_vertex_property("string")
    g.vertex_properties["internal_ids"] = v_internal_id

    v_id = g.new_vertex_property("string")
    g.vertex_properties["id"] = v_id

    # For protein nodes
    v_reviewed = g.new_vertex_property("boolean")
    g.vertex_properties["reviewed"] = v_reviewed

    # store nodes to connect them when creating edges
    vertices = {}
    drug_vertices = {}
    # add vertices

    print(f'loading nodes for {identifier}')

    is_entrez = (identifier == 'entrez' or identifier == 'ncbigene')
    is_symbol = identifier == 'symbol'
    is_uniprot = identifier == 'uniprot'
    is_ensg = (identifier == 'ensg' or identifier == 'ensembl')
    is_internal = identifier == None

    if is_ensg or is_internal:
        ensembl_set = defaultdict(set)
        for node in models.EnsemblGene.objects.all():
            ensembl_set[node.protein_id].add(node.name)

    node_id_map = defaultdict(set)
    drugstone_id_to_node = dict()

    if isReviewed:
        proteins = models.Protein.objects.filter(isReviewed=True)
    else:
        proteins = models.Protein.objects.all()

    for node in proteins:
        if is_entrez:
            if len(node.entrez) != 0:
                node_id_map[node.entrez].add(node.id)
                drugstone_id_to_node[node.id].add(node)
        elif is_symbol:
            if len(node.gene) != 0:
                node_id_map[node.gene].add(node.id)
                drugstone_id_to_node[node.id] = node
        elif is_uniprot:
            node_id_map[node.uniprot_code].add(node.id)
            drugstone_id_to_node[node.id] = node
        elif is_ensg:
            for id in ensembl_set[node.id]:
                node_id_map[id].add(node.id)
                drugstone_id_to_node[node.id] = node
        elif is_internal:
            node_id_map[node.id].add(node.id)
            drugstone_id_to_node[node.id] =node

            v_uniprot = g.new_vertex_property("string")
            g.vertex_properties["uniprot"] = v_uniprot

            v_symbol = g.new_vertex_property("string")
            g.vertex_properties["symbol"] = v_symbol

            v_entrez = g.new_vertex_property("string")
            g.vertex_properties["entrez"] = v_entrez

            v_ensembl = g.new_vertex_property("string")
            g.vertex_properties["ensembl"] = v_ensembl

        for id, internal_ids in node_id_map.items():
            v = g.add_vertex()
            v_type[v] = 'protein'
            v_internal_id[v] = ",".join({f"p{id}" for id in internal_ids})
            for drugstone_id in internal_ids:
                node = drugstone_id_to_node[drugstone_id]
                v_reviewed[v] = node.isReviewed
                v_name[v] = node.gene
                vertices[drugstone_id] = v
            if is_internal:
                node = drugstone_id_to_node[id]
                v_uniprot[v] = node.uniprot_code
                v_symbol[v] = node.gene
                v_entrez[v] = node.entrez
                if node.id in ensembl_set.keys():
                    v_ensembl[v] = ",".join({f"{id}" for id in ensembl_set[node.id]})


    uniq_edges = set()

    n = 0
    for edge_raw in _internal_ppis(dataset):
        id1 = edge_raw.from_protein_id
        id2 = edge_raw.to_protein_id
        if id1 > id2:
            tmp = id1
            id1 = id2
            id2 = tmp
        hash = f'{id1}_{id2}'
        if hash not in uniq_edges and id1 in vertices and id2 in vertices:
            uniq_edges.add(hash)
            e = g.add_edge(vertices[id1], vertices[id2])
            n += 1
            e_type[e] = 'protein-protein'
            e_directed[e] = edge_raw.is_directed
            e_inhibition[e] = edge_raw.is_inhibition
            e_stimulation[e] = edge_raw.is_stimulation
    print("done with PPI edges: ", n)

    # remove unconnected proteins
    delete_vertices = set()
    for vertex in vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)
    print("removing unconnected proteins: ", len(delete_vertices))

    # remove unconnected drugs
    for vertex in drug_vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)

    g.remove_vertex(reversed(sorted(delete_vertices)), fast=True)
    Path('./data/Networks/').mkdir(parents=True, exist_ok=True)
    g.save(filename, fmt=fmt)
    print(f"Created file {filepath}")
    print("Size of graph - nodes: ", g.num_vertices(), " edges: ", g.num_edges())
    return filepath

def create_gt(params: List[str]) -> None:
    """Fetches all required information to build a graph-tools file for given
    PPI and PDI dataset names (params). Builds the graph-tools file and saves it in 
    the data/Networks folder.

    Args:
        params (Tuple[str, str]): Protein-protein-dataset name, Protein-drug-dataset name
    """
    ppi_dataset, pdi_dataset, identifier, isReviewed = params

    licensed = ppi_dataset.licenced or pdi_dataset.licenced
    # get data from api
    
        # save graph
    filename = f"./data/Networks/{identifier}_{ppi_dataset.name}-{pdi_dataset.name}"
    if licensed:
        filename += "_licenced"

    if isReviewed:
        filename += "_reviewed"

    filename += ".gt"
    
    print(f'Creating {filename}')

    g = gt.Graph(directed=False)

    e_type = g.new_edge_property("string")

    v_type = g.new_vertex_property("string")
    v_name = g.new_vertex_property("string")

    # for drugs
    v_status = g.new_vertex_property("string")
    v_drug_id = g.new_vertex_property("string")
    v_internal_id = g.new_vertex_property("string")


    g.edge_properties["type"] = e_type
    # g.edge_properties["drugstone_id"] = e_type

    g.vertex_properties["type"] = v_type
    g.vertex_properties["name"] = v_name
    g.vertex_properties["status"] = v_status
    g.vertex_properties["drug_id"] = v_drug_id
    g.vertex_properties["internal_id"] = v_internal_id

    # store nodes to connect them when creating edges
    vertices = {}
    drug_vertices = {}
    # add vertices

    print(f'loading nodes for {identifier}')

    is_entrez = (identifier == 'entrez' or identifier == 'ncbigene')
    is_symbol = identifier == 'symbol'
    is_uniprot = identifier == 'uniprot'
    is_ensg = (identifier == 'ensg' or identifier == 'ensembl')

    if is_ensg:
        ensembl_set = defaultdict(set)
        for node in models.EnsemblGene.objects.all():
            ensembl_set[node.protein_id].add(node.name)

    node_id_map = defaultdict(set)
    drugstone_ids_to_node_ids = defaultdict(set)

    if isReviewed:
        proteins = models.Protein.objects.filter(isReviewed=True)
    else:
        proteins = models.Protein.objects.all()

    for node in proteins:
        if is_entrez:
            if len(node.entrez) != 0:
                node_id_map[node.entrez].add(node.id)
                drugstone_ids_to_node_ids[node.id].add(node.entrez)
        elif is_symbol:
            if len(node.gene) != 0:
                node_id_map[node.gene].add(node.id)
                drugstone_ids_to_node_ids[node.id].add(node.gene)
        elif is_uniprot:
            node_id_map[node.uniprot_code].add(node.id)
            drugstone_ids_to_node_ids[node.id].add(node.uniprot_code)
        elif is_ensg:
            for id in ensembl_set[node.id]:
                node_id_map[id].add(node.id)
                drugstone_ids_to_node_ids[node.id].add(id)

    for id, nodes in node_id_map.items():
        v = g.add_vertex()
        v_type[v] = 'protein'
        v_internal_id[v] = id
        for drugstone_id in nodes:
            vertices[drugstone_id] = v

    for node in models.Drug.objects.all():
        v = g.add_vertex()
        v_type[v] = 'drug'
        v_status[v] = node.status
        v_internal_id[v] = f'dr{node.id}'

        drug_vertices[node.id] = v


    uniq_edges = set()

    n = 0
    for edge_raw in _internal_ppis(ppi_dataset):
        id1 = edge_raw.from_protein_id
        id2 = edge_raw.to_protein_id
        if id1 > id2:
            tmp = id1
            id1 = id2
            id2 = tmp
        hash = f'{id1}_{id2}'
        if hash not in uniq_edges and id1 in vertices and id2 in vertices:
            uniq_edges.add(hash)
            e = g.add_edge(vertices[id1], vertices[id2])
            n += 1
            e_type[e] = 'protein-protein'
    print("done with PPI edges: ", n)

    uniq_edges = set()
    n = 0
    print(f'loading drug_edges/{pdi_dataset}')
    for edge_raw in _internal_pdi(pdi_dataset):
        id1 = edge_raw.drug_id
        id2 = edge_raw.protein_id
        hash = f'{id1}_{id2}'
        if hash not in uniq_edges and id1 in drug_vertices and id2 in vertices:
            uniq_edges.add(hash)
            e = g.add_edge(drug_vertices[id1], vertices[id2])
            n += 1
            e_type[e] = 'drug-protein'
    print("done with drug edges: ", n)

    # remove unconnected proteins
    delete_vertices = set()
    for vertex in vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)
    print("removing unconnected proteins: ", len(delete_vertices))

    # remove unconnected drugs
    for vertex in drug_vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)

    g.remove_vertex(reversed(sorted(delete_vertices)), fast=True)
    Path('./data/Networks/').mkdir(parents=True, exist_ok=True)
    g.save(filename)
    print(f"Created file {filename}")
    print("Size of graph - nodes: ", g.num_vertices(), " edges: ", g.num_edges())
    return


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **kwargs):
        ppi_datasets = models.PPIDataset.objects.all()

        pdi_datasets = models.PDIDataset.objects.all()

        licenced_ppi_dataset = {ppi.name: ppi for ppi in ppi_datasets if ppi.licenced}
        licenced_pdi_dataset = {pdi.name: pdi for pdi in pdi_datasets if pdi.licenced}

        uniq_combis = set()
        parameter_combinations = []
        for protein_interaction_dataset in ppi_datasets:
            for pdi_dataset in pdi_datasets:
                ppi_ds = protein_interaction_dataset
                pdi_ds = pdi_dataset
                licenced = ppi_ds.licenced or pdi_ds.licenced
                if licenced:
                    ppi_ds = licenced_ppi_dataset[
                        ppi_ds.name] if protein_interaction_dataset.name in licenced_ppi_dataset else ppi_ds
                    pdi_ds = licenced_pdi_dataset[
                        pdi_ds.name] if pdi_ds.name in licenced_pdi_dataset else pdi_ds
                hash = f'{ppi_ds.name}-{pdi_ds.name}_{licenced}'
                if hash in uniq_combis:
                    continue
                uniq_combis.add(hash)
                for identifier in ['ensg', 'symbol', 'entrez', 'uniprot']:
                    for isReviewed in [True, False]:
                        parameter_combinations.append([ppi_ds, pdi_ds, identifier, isReviewed])
        # close all database connections so subprocesses will create their own connections
        # this prevents the processes from running into problems because of using the same connection
        db.connections.close_all()
        pool = multiprocessing.Pool(KERNEL)
        pool.map(create_gt, parameter_combinations)
