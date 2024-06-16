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


def create_gt(params: List[str]) -> None:
    """Fetches all required information to build a graph-tools file for given
    PPI and PDI dataset names (params). Builds the graph-tools file and saves it in 
    the data/Networks folder.

    Args:
        params (Tuple[str, str]): Protein-protein-dataset name, Protein-drug-dataset name
    """
    ppi_dataset, pdi_dataset, identifier = params

    licensed = ppi_dataset.licenced or pdi_dataset.licenced
    # get data from api
    
        # save graph
    filename = f"./data/Networks/{identifier}_{ppi_dataset.name}-{pdi_dataset.name}"
    if licensed:
        filename += "_licenced"
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

    for node in models.Protein.objects.all():
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
    print("done with nodes")

    print(f"adding drugs")
    for node in models.Drug.objects.all():
        v = g.add_vertex()
        v_type[v] = 'drug'
        v_status[v] = node.status
        v_internal_id[v] = f'dr{node.id}'

        drug_vertices[node.id] = v

    print("done with drugs")

    # add edges
    print(f'adding ppi_edges/{ppi_dataset}')

    uniq_edges = set()

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
            e_type[e] = 'protein-protein'
    print("done with edges")

    uniq_edges = set()

    print(f'loading drug_edges/{pdi_dataset}')
    for edge_raw in _internal_pdis(pdi_dataset):
        id1 = edge_raw.drug_id
        id2 = edge_raw.protein_id
        hash = f'{id1}_{id2}'
        if hash not in uniq_edges and id1 in drug_vertices and id2 in vertices:
            uniq_edges.add(hash)
            e = g.add_edge(drug_vertices[id1], vertices[id2])
            e_type[e] = 'drug-protein'
    print("done with drug edges")

    # remove unconnected proteins
    delete_vertices = set()
    for vertex in vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)

    # remove unconnected drugs
    for vertex in drug_vertices.values():
        if vertex.out_degree() == 0:
            delete_vertices.add(vertex)

    g.remove_vertex(reversed(sorted(delete_vertices)), fast=True)
    Path('./data/Networks/').mkdir(parents=True, exist_ok=True)
    g.save(filename)
    print(f"Created file {filename}")
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
                    parameter_combinations.append([ppi_ds, pdi_ds, identifier])
        # close all database connections so subprocesses will create their own connections
        # this prevents the processes from running into problems because of using the same connection
        db.connections.close_all()
        pool = multiprocessing.Pool(KERNEL)
        pool.map(create_gt, parameter_combinations)
