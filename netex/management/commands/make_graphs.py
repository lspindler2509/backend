from typing import List, Tuple
import graph_tool.all as gt
from netex import models
from netex import serializers
import json
import multiprocessing
from django import db

from django.core.management import BaseCommand
import django

django.setup()

KERNEL = 6

def _internal_expression_scores(netex_id: str) -> dict:
    """ Looks up the tissue specific expression scores for a given protein.
    The scores are loaded from the django database.

    Args:
        netex_id (str): netex id of protein in format 'pxxxx'

    Returns:
        dict: keys are tissue-names and values are the expression scores
    """
    protein_object = models.Protein.objects.get(id=int(netex_id[1:]))

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

def _internal_pdis(dataset_name: str) -> List[dict]:
    """ Fetches all internal protein-drug interactions for a given dataset.
    Interactions are taken from the django database.

    Args:
        dataset_name (str): Name of the dataset, e.g. "DrugBank"

    Returns:
        List[dict]: List of representaions of interaction objects
    """
    # get all interactions
    node_node_interaction_objects = models.ProteinDrugInteraction.objects.filter(
        pdi_dataset__name=dataset_name
    )
    node_node_interactions = serializers.ProteinDrugInteractionSerializer(many=True) \
        .to_representation(node_node_interaction_objects)

    return node_node_interactions

def _internal_ppis(dataset_name: str) -> List[dict]:
    """ Fetches all internal protein-protein interactions for a given dataset.
    Interactions are taken from the django database.

    Args:
        dataset_name (str): Name of the dataset, e.g. "BioGRID"

    Returns:
        List[dict]: List of representaions of interaction objects
    """
    # get all interactions
    node_node_interaction_objects = models.ProteinProteinInteraction.objects.filter(
        ppi_dataset__name=dataset_name
    )
    node_node_interactions = serializers.ProteinProteinInteractionSerializer(many=True) \
        .to_representation(node_node_interaction_objects)

    return node_node_interactions


def create_gt(params: Tuple[str, str]) -> None:
    """Fetches all required information to build a graph-tools file for given
    PPI and PDI dataset names (params). Builds the graph-tools file and saves it in 
    the data-NetExpander/Networks folder.

    Args:
        params (Tuple[str, str]): Protein-protein-dataset name, Protein-drug-dataset name
    """
    ppi_dataset, pdi_dataset = params
    # get data from api
    data = {}

    print(f'loading nodes')
    data['nodes'] = serializers.ProteinSerializer(many=True).to_representation(
        models.Protein.objects.all()
    ) 

    print(f'loading edges/{ppi_dataset}')
    data['edges'] = _internal_ppis(ppi_dataset)

    print(f'loading drugs')
    data['drugs'] = serializers.DrugSerializer(many=True).to_representation(
        models.Drug.objects.all()
    ) 
    print(f'loading drug_edges/{pdi_dataset}')
    data['drug_edges'] = _internal_pdis(pdi_dataset)

    g = gt.Graph(directed=False)
    e_type = g.new_edge_property("string")

    v_type = g.new_vertex_property("string")
    v_name = g.new_vertex_property("string")
    v_netex_id = g.new_vertex_property("string")
    v_entrez = g.new_vertex_property("string")
    v_expression = g.new_vertex_property("string")

    # for drugs
    v_status = g.new_vertex_property("string")
    v_drug_id = g.new_vertex_property("string")

    g.edge_properties["type"] = e_type
    g.edge_properties["netex_id"] = e_type

    g.vertex_properties["type"] = v_type
    g.vertex_properties["name"] = v_name
    g.vertex_properties["netex_id"] = v_netex_id
    g.vertex_properties["entrez"] = v_entrez
    g.vertex_properties["status"] = v_status
    g.vertex_properties["drug_id"] = v_drug_id
    g.vertex_properties["expression"] = v_expression

    # store nodes to connect them when creating edges
    vertices = {}
    drug_vertices = {}
    # add vertices
    print("adding nodes")
    # extend node data by cancer nodes, we create a normal node for each cancer node.
    # on reading the data, we decide which one to keep based on the user selected cancer types
    for node in data['nodes']:
        v = g.add_vertex()
        v_type[v] = 'protein'
        v_name[v] = node['symbol']
        v_netex_id[v] = node['netex_id']
        v_entrez[v] = node['entrez']
        v_expression[v] = json.dumps(_internal_expression_scores(node['netex_id']))

        vertices[node['netex_id']] = v

    print("done with nodes")

    print(f"adding ${len(data['drugs'])} drugs")
    for node in data['drugs']:
        v = g.add_vertex()
        v_type[v] = 'drug'
        v_name[v] = node['label']
        v_status[v] = node['status']
        v_netex_id[v] = node['netex_id']
        v_drug_id[v] = node['drug_id']

        drug_vertices[node['netex_id']] = v
    print("done with drugs")

    # add edges
    print("adding edges")
    for edge in data['edges']:
        a = vertices[edge['protein_a']]
        b = vertices[edge['protein_b']]
        e = g.add_edge(a, b)
        e_type[e] = 'protein-protein'
    print("done with edges")

    print("adding drug edges")
    for edge in data['drug_edges']:
        protein, drug = vertices[edge['protein']], drug_vertices[edge['drug']]
        e = g.add_edge(drug, protein)
        e_type[e] = 'drug-protein'
    print("done with drug edges")

    # save graph
    filename = f"./data-NetExpander/Networks/internal_{ppi_dataset}_{pdi_dataset}.gt"
    g.save(filename)
    print(f"Created file {filename}")
    return


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **kwargs):
        ppi_datasets = models.PPIDataset.objects.all()
        ppi_datasets_names = [e.name for e in ppi_datasets]

        pdi_datasets = models.PDIDataset.objects.all()
        pdi_datasets_names = [e.name for e in pdi_datasets]

        parameter_combinations = []
        for protein_interaction_dataset in ppi_datasets_names:
            for pdi_dataset in pdi_datasets_names:
                parameter_combinations.append((protein_interaction_dataset, pdi_dataset))

        # close all database connections so subprocesses will create their own connections
        # this prevents the processes from running into problems because of using the same connection
        db.connections.close_all()
        pool = multiprocessing.Pool(KERNEL)
        pool.map(create_gt, parameter_combinations)
