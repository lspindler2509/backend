from typing import List, Tuple
import graph_tool.all as gt
from drugstone import models
import multiprocessing
from django import db

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


def _internal_pdis(dataset) -> List[models.ProteinDrugInteraction]:
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


def create_gt(params: Tuple) -> None:
    """Fetches all required information to build a graph-tools file for given
    PPI and PDI dataset names (params). Builds the graph-tools file and saves it in 
    the data/Networks folder.

    Args:
        params (Tuple[str, str]): Protein-protein-dataset name, Protein-drug-dataset name
    """
    ppi_dataset, pdi_dataset = params
    licensed = ppi_dataset.licenced or pdi_dataset.licenced
    # get data from api

    g = gt.Graph(directed=False)
    e_type = g.new_edge_property("string")

    v_type = g.new_vertex_property("string")
    v_name = g.new_vertex_property("string")
    v_drugstone_id = g.new_vertex_property("string")
    v_entrez = g.new_vertex_property("string")
    v_expression = g.new_vertex_property("string")

    # for drugs
    v_status = g.new_vertex_property("string")
    v_drug_id = g.new_vertex_property("string")

    g.edge_properties["type"] = e_type
    g.edge_properties["drugstone_id"] = e_type

    g.vertex_properties["type"] = v_type
    g.vertex_properties["name"] = v_name
    g.vertex_properties["drugstone_id"] = v_drugstone_id
    g.vertex_properties["entrez"] = v_entrez
    g.vertex_properties["status"] = v_status
    g.vertex_properties["drug_id"] = v_drug_id
    g.vertex_properties["expression"] = v_expression

    # store nodes to connect them when creating edges
    vertices = {}
    drug_vertices = {}
    # add vertices

    # print("adding nodes")
    print(f'loading nodes')
    # extend node data by cancer nodes, we create a normal node for each cancer node.
    # on reading the data, we decide which one to keep based on the user selected cancer types
    for node in models.Protein.objects.all():
        v = g.add_vertex()
        v_type[v] = 'protein'
        v_drugstone_id[v] = f"p{node.id}"

        vertices[node.id] = v

    print("done with nodes")

    print(f"adding drugs")
    for node in models.Drug.objects.all():
        v = g.add_vertex()
        v_type[v] = 'drug'
        v_status[v] = node.status
        v_drugstone_id[v] = f'dr{node.id}'

        drug_vertices[node.id] = v
    print("done with drugs")

    # add edges
    print(f'adding ppi_edges/{ppi_dataset}')
    for edge_raw in _internal_ppis(ppi_dataset):
        e = g.add_edge(vertices[edge_raw.from_protein_id], vertices[edge_raw.to_protein_id])
        e_type[e] = 'protein-protein'
    print("done with edges")


    print(f'loading drug_edges/{pdi_dataset}')
    for edge_raw in _internal_pdis(pdi_dataset):
        e = g.add_edge(drug_vertices[edge_raw.drug_id], vertices[edge_raw.protein_id])
        e_type[e] = 'drug-protein'
    print("done with drug edges")

    # remove unconnected proteins
    delete_vertices = set()
    for vertex in vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)

    #remove unconnected drugs
    for vertex in drug_vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)

    g.remove_vertex(reversed(sorted(delete_vertices)), fast=True)

    # save graph
    filename = f"./data/Networks/internal_{ppi_dataset.name}_{pdi_dataset.name}"
    if licensed:
        filename += "_licenced"
    filename += ".gt"
    g.save(filename)
    print(f"Created file {filename}")
    return


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **kwargs):
        ppi_datasets = models.PPIDataset.objects.all()

        pdi_datasets = models.PDIDataset.objects.all()

        parameter_combinations = []
        for protein_interaction_dataset in ppi_datasets:
            for pdi_dataset in pdi_datasets:
                parameter_combinations.append((protein_interaction_dataset, pdi_dataset))

        # close all database connections so subprocesses will create their own connections
        # this prevents the processes from running into problems because of using the same connection
        db.connections.close_all()
        pool = multiprocessing.Pool(KERNEL)
        pool.map(create_gt, parameter_combinations)
