import csv
import math
import random
import string
import time
import uuid
import mimetypes
from collections import defaultdict
from typing import List

import pandas as pd
import networkx as nx
from django.http import HttpResponse, JsonResponse
from django.db.models import Q, Max
from django.db import IntegrityError
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import parsers, views
from rest_framework.views import APIView
import graph_tool as gt
import networkx as nx
from rest_framework.response import Response
from django.utils.encoding import smart_str
from django.http import StreamingHttpResponse
from wsgiref.util import FileWrapper

from drugstone.util.mailer import bugreport
from drugstone.util.property_calulations import calculate_properties
from drugstone.util.query_db import (
    fetch_edges_from_input,
    map_edges,
    query_proteins_by_identifier,
    clean_proteins_from_compact_notation,
    fetch_node_information,
    update_result,
)

from drugstone.models import *
from drugstone.serializers import *
from drugstone.backend_tasks import (
    start_task,
    refresh_from_redis,
    task_stats,
    task_result,
    task_parameters,
)

from tasks.pathway_enrichment import get_all_node_scores, parse_pathway;
from tasks.create_genesets import parse_genesets;

from drugstone.settings import DEFAULTS
import os
from tasks.util.custom_network import remove_ppi_edges


def get_ppi_ds(source, licenced):
    ds = models.PPIDataset.objects.filter(name__iexact=source, licenced=licenced).last()
    if ds is None and licenced:
        return get_ppi_ds(source, False)
    return ds


def get_pdi_ds(source, licenced):
    ds = models.PDIDataset.objects.filter(name__iexact=source, licenced=licenced).last()
    if ds is None and licenced:
        return get_pdi_ds(source, False)
    return ds


def get_pdis_ds(source, licenced):
    ds = models.PDisDataset.objects.filter(
        name__iexact=source, licenced=licenced
    ).last()
    if ds is None and licenced:
        return get_pdis_ds(source, False)
    return ds


def get_drdis_ds(source, licenced):
    ds = models.DrDiDataset.objects.filter(
        name__iexact=source, licenced=licenced
    ).last()
    if ds is None and licenced:
        return get_drdis_ds(source, False)
    return ds


class FileUploadView(views.APIView):
    parser_classes = [parsers.MultiPartParser]

    def post(self, request, filename, format=None):
        file_obj = request.data['file']
        try:
            parsed_network = self.parseFile(file_obj)
            return Response(parsed_network)
        except Exception as e:
            print("Error occured during file parsing: ", e)
            return Response(False)

    def parseFile(self, file):
        if file.name.endswith('.graphml'):
            file_content = file.read().decode('utf-8')
            G = nx.parse_graphml(file_content)
            nodes = []
            for node in G.nodes():
                node_id = G.nodes[node].get("name", str(node))
                group_value = G.nodes[node].get('group', 'default')

                node_data = {'id': node_id, 'group': group_value}
                node_data['properties'] = {key: value for key, value in G.nodes[node].items()}

                nodes.append(node_data)

            edges = [{'from': str(edge[0]), 'to': str(edge[1])} for edge in G.edges()]
            return {'nodes': nodes, 'edges': edges}

        if file.name.endswith('.gt'):
            g = gt.load_graph(file, fmt="gt")
            nodes = []
            hasGroup: bool = "group" in g.vertex_properties
            for node in g.vertices():
                node_data = {'id': str(g.vertex_properties["name"][node])}
                node_data["properties"] = {}
                for prop_name, prop_map in g.vertex_properties.items():
                    if prop_name != "name":
                        node_data["properties"][prop_name] = prop_map[node]
                node_data["group"] = g.vertex_properties["group"][node] if hasGroup else "default"
                nodes.append(node_data)
            edges = [
                {'from': g.vertex_properties["name"][edge.source()], 'to': g.vertex_properties["name"][edge.target()]}
                for edge in g.edges()]
            return {'nodes': nodes, 'edges': edges}

        nodes = []
        edges = []
        unique_nodes = set()
        unique_edges = set()

        contents = file.read().decode('utf-8').splitlines()

        for line in contents:
            line = line.strip()

            if not line:
                continue

            if file.name.endswith('.csv'):
                clean_from, clean_to = line.split(',')
                clean_from = clean_from.strip().split('.')[1] if len(
                    clean_from.strip().split('.')) > 1 else clean_from.strip()
                clean_to = clean_to.strip().split('.')[1] if len(clean_to.strip().split('.')) > 1 else clean_to.strip()

            elif file.name.endswith('.sif'):
                parts = line.split()
                if len(parts) != 3:
                    isolated_node_id = parts[0].strip().split('.')[1] if len(parts[0].strip().split('.')) > 1 else \
                        parts[0].strip()
                    if not isolated_node_id in unique_nodes:
                        nodes.append({'id': isolated_node_id, 'group': "default"})
                    unique_nodes.add(isolated_node_id)
                    continue

                clean_from, _, clean_to = parts
                clean_from = clean_from.strip().split('.')[1] if len(
                    clean_from.strip().split('.')) > 1 else clean_from.strip()
                clean_to = clean_to.strip().split('.')[1] if len(clean_to.strip().split('.')) > 1 else clean_to.strip()

            else:
                return Response(False)

            edge_key = '-'.join(sorted([clean_from, clean_to]))

            if clean_from not in unique_nodes:
                nodes.append({'id': clean_from, 'group': "default"})
                unique_nodes.add(clean_from)
            if clean_to not in unique_nodes:
                nodes.append({'id': clean_to, 'group': "default"})
                unique_nodes.add(clean_to)

            if edge_key not in unique_edges:
                edges.append({'from': clean_from, 'to': clean_to})
                unique_edges.add(edge_key)

        return {'nodes': nodes, 'edges': edges}


class TaskView(APIView):
    def post(self, request) -> Response:
        chars = string.ascii_lowercase + string.ascii_uppercase + string.digits
        token_str = "".join(random.choice(chars) for _ in range(32))
        parameters = request.data["parameters"]
        licenced = parameters.get("licenced", False)
        algorithm = request.data["algorithm"]

        # find databases based on parameter strings
        parameters["ppi_dataset"] = PPIDatasetSerializer().to_representation(
            get_ppi_ds(parameters.get("ppi_dataset", DEFAULTS["ppi"]), licenced)
        )

        parameters["pdi_dataset"] = PDIDatasetSerializer().to_representation(
            get_pdi_ds(parameters.get("pdi_dataset", DEFAULTS["pdi"]), licenced)
        )

        # if algorithm in ['connect', 'connectSelected', 'quick', 'super']:
        #     parameters["num_trees"] = 5
        #     parameters["tolerance"] = 5
        #     parameters["hub_penalty"] = 0.5

        task = Task.objects.create(
            token=token_str,
            target=request.data["target"],
            algorithm=algorithm,
            parameters=json.dumps(parameters),
        )
        start_task(task)
        task.save()

        return Response(
            {
                "token": token_str,
            }
        )

    def get(self, request) -> Response:
        token_str = request.query_params["token"]
        task = Task.objects.get(token=token_str)

        if not task.done and not task.failed:
            refresh_from_redis(task)
            task.save()

        return Response(
            {
                "token": task.token,
                "info": TaskSerializer().to_representation(task),
                "stats": task_stats(task),
            }
        )


@api_view(["GET"])
def get_license(request) -> Response:
    from drugstone.management.includes.DatasetLoader import import_license
    return Response({"license": import_license()})


@api_view(["POST"])
def create_genesets(request) -> Response:
    kegg = request.query_params["kegg"]
    reactome = request.query_params["reactome"]
    wiki = request.query_params["wiki"]
    print("Creating genesets")
    parse_genesets(kegg, reactome, wiki, False)
    print("Created genesets unreviewed")
    parse_genesets(kegg, reactome, wiki, True)
    print("Created genesets reviewed")
    return Response("worked!")


@api_view(["GET"])
def get_default_params(request) -> Response:
    algorithm = request.GET.get("algorithm")
    connect = {
        "algorithm": "multisteiner",
        "numTrees": 5,
        "tolerance": 5,
        "hubPenalty": 0.5,
    }
    quick = {
        "algorithm": "closeness",
        "result_size": 50,
        "hub_penalty": 0,
        "include_non_approved_drugs": False,
        "include_indirect_drugs": False,
    }
    resp = {}
    if algorithm in ["quick", "super", "connect", "connectSelected"]:
        resp["protein"] = connect
    if algorithm in ["quick", "super"]:
        resp["drug"] = quick
    return Response(resp)


@api_view(["POST"])
def fetch_edges(request) -> Response:
    """Retrieves interactions between nodes given as a list of drugstone IDs.

    Args:
        request (HttpRequest): With keys 'nodes' representing nodes and 'dataset' representing the
        protein-protein interaction dataset.

    Returns:
        Response: List of edges which are objects with 'from' and to ' attribtues'
    """
    dataset = request.data.get("dataset", DEFAULTS["ppi"])
    drugstone_ids = set()
    for node in request.data.get("nodes", "[]"):
        if "drugstone_id" in node:
            if isinstance(node["drugstone_id"], list):
                for id in node["drugstone_id"]:
                    drugstone_ids.add(id[1:])
            else:
                drugstone_ids.add(node["drugstone_id"])
    licenced = request.data.get("licenced", False)
    dataset_object = get_ppi_ds(dataset, licenced)
    interaction_objects = models.ProteinProteinInteraction.objects.filter(
        Q(ppi_dataset=dataset_object)
        & Q(from_protein__in=drugstone_ids)
        & Q(to_protein__in=drugstone_ids)
    )

    return Response(
        ProteinProteinInteractionSerializer(many=True).to_representation(
            interaction_objects
        )
    )


@api_view(['GET'])
def searchProteins(request) -> Response:
    try:
        query = request.query_params.get("query", "")
        limit = request.query_params.get("limit", 20)
        identifier = request.query_params.get("identifier", "symbol")
        label = request.query_params.get("label", "")
        reviewed = request.query_params.get("reviewed", False)
        reviewed = True if reviewed == "true" else False
        if not query:
            return Response([])

        # Filter proteins by uniprot_code, gene (symbol), entrez, or related EnsemblGene name
        proteins = models.Protein.objects.filter(
            Q(uniprot_code__icontains=query) |
            Q(gene__icontains=query) |
            Q(entrez__icontains=query) |
            Q(ensg__name__icontains=query)
        )

        if reviewed:
            proteins = proteins.filter(isReviewed=True)
        proteins = proteins.distinct()[:limit]

        uniprot_ids = list(proteins.values_list("uniprot_code", flat=True))
        mapped_nodes, _ = query_proteins_by_identifier(uniprot_ids, "uniprot", reviewed)

        for node in mapped_nodes:
            node["label"] = node[label][0] if label in node and node[label] else node["uniprot"][0]
            node["id"] = node[identifier][0] if identifier in node and node[identifier] else node["uniprot"][0]

        return Response(mapped_nodes)
    except Exception as e:
        print("An error occured while searching for proteins: ", e)
        return Response([])


@api_view(["POST"])
def convert_compact_ids(request) -> Response:
    nodes = request.data.get("nodes", "[]")
    identifier = request.data.get("identifier", "")
    cleaned = clean_proteins_from_compact_notation(nodes, identifier)
    return Response(cleaned)


@api_view(["POST"])
def prepare_pruning(request) -> Response:
    try:
        data = json.loads(request.body)
    except json.JSONDecodeError:
        return JsonResponse({"error": "Invalid JSON"}, status=400)
    try:
        nodes = data.get("nodes", [])
        pruning_attribute = data.get("pruning_attribute", "")
        pruning_result = {}

        for node in nodes:
            properties = node.get("properties", {})
            value = properties.get(pruning_attribute)
            if isinstance(value, str):
                pruning_result.setdefault("unique_values", set()).add(value)
                pruning_result["type"] = "string"
            elif isinstance(value, (int, float)):
                pruning_result["min"] = math.floor(min(pruning_result.get("min", value), value))
                pruning_result["max"] = math.ceil(max(pruning_result.get("max", value), value))
                if not pruning_result.get("type", False) or pruning_result["type"] == "int":
                    pruning_result["type"] = type(value).__name__

        if "unique_values" in pruning_result:
            pruning_result["unique_values"] = list(pruning_result["unique_values"])

        return Response(pruning_result)
    except Exception as e:
        print("An error occured while preparing pruning: ", e)
        return Response({})


@api_view(["POST"])
def recalculate_statistics(request) -> Response:
    try:
        data = json.loads(request.body)
        network = data.get("network", {})
        nodes = network.get("nodes", [])
        edges = network.get("edges", [])
        config = data.get("config", {})
    except json.JSONDecodeError as e:
        print("Something went wrong while parsing the body!", e)
        return JsonResponse({"error": "Invalid JSON"}, status=400)

    calculateProperties = config.get("calculateProperties", False)
    if not calculateProperties:
        return Response(calculate_properties(nodes, None, None, None, False))
    id_space = config.get("identifier", "symbol")
    custom_edges = config.get("custom_edges", False)
    no_default_edges = config.get("exclude_drugstone_ppi_edges", False)
    ppi_dataset = config.get("interactionProteinProtein")
    pdi_dataset = config.get("interactionDrugProtein")

    filename = f"{id_space}_{ppi_dataset}-{pdi_dataset}"
    if config.get("licensedDatasets", False):
        filename += "_licenced"
    if config.get("reviewed", False):
        filename += "_reviewed"
    filename = os.path.join("./data/Networks/", filename + ".gt")
    graph = gt.load_graph(filename)
    if custom_edges:
        if no_default_edges:
            # clear all edges with type "protein-protein"
            graph = remove_ppi_edges(graph)
        edges = edges
        graph = add_edges(graph, edges)

    return Response(calculate_properties(nodes, graph, id_space, edges))


@api_view(["POST"])
def overlay_directed_edges(request) -> Response:
    try:
        data = json.loads(request.body)
        ppi_dataset = data.get("ppi_dataset", "")
        licenced = data.get("licenced", False)
        ppi_dataset = PPIDatasetSerializer().to_representation(get_ppi_ds(ppi_dataset, licenced))
        edges = data.get("edges", [])
        nodes_mapped_dict = data.get("nodes_mapped_dict", {})
        drugstone_mapping = data.get("drugstone_mapping", False)
    except json.JSONDecodeError:
        return JsonResponse({"error": "Invalid JSON"}, status=400)

    edges_overlayed = map_edges(ppi_dataset, edges, nodes_mapped_dict, drugstone_mapping, "drugstoneId")
    edges_with_ids = []
    edge_id_map = {(edge["from"], edge["to"]): edge["id"] for edge in edges}
    edge_id_map.update({(edge["to"], edge["from"]): edge["id"] for edge in edges})
    for edge in edges_overlayed:
        edge.pop("groupName", None)
        edge_id = edge_id_map.get((edge["from"], edge["to"]))
        if edge_id:
            edge["id"] = edge_id
        edges_with_ids.append(edge)

    return Response(edges_with_ids)


@api_view(["POST"])
def prune(request) -> Response:
    try:
        data = json.loads(request.body)
        network = data.get("network", {})
        nodes = network.get("nodes", [])
        edges = network.get("edges", [])
        pruning_attribute = data.get("pruning_attribute", "")
        prune_orphan_nodes = data.get("pruneOrphanNodes", False)
        cutoff = data.get("cutoff", None)
        pruningDirection = data.get("pruningDirection", "greater")
        unique_values = data.get("unique_values", [])
    except json.JSONDecodeError:
        return JsonResponse({"error": "Invalid JSON"}, status=400)

    pruned_node_ids = set()

    if len(unique_values) > 0:
        unique_values_set = {value for value in unique_values if value}
        pruned_node_ids = {node["id"] for node in nodes if
                           node["properties"].get(pruning_attribute, "") in unique_values_set}
    elif cutoff is not None:
        if pruningDirection == "greater":
            pruned_node_ids = {node["id"] for node in nodes if
                               (node["properties"].get(pruning_attribute, cutoff - 1) >= cutoff)}
        elif pruningDirection == "lesser":
            pruned_node_ids = {node["id"] for node in nodes if
                               (node["properties"].get(pruning_attribute, cutoff + 1) <= cutoff)}

    pruned_edges = [
        edge for edge in edges
        if edge.get("from") in pruned_node_ids and edge.get("to") in pruned_node_ids
    ]

    if prune_orphan_nodes:
        connected_node_ids = {edge["from"] for edge in pruned_edges} | {edge["to"] for edge in pruned_edges}
        orphan_node_ids = {node["id"] for node in nodes if node["id"] not in connected_node_ids}
        pruned_node_ids = pruned_node_ids - orphan_node_ids

    for node in nodes:
        if node["id"] not in pruned_node_ids:
            node["to_be_pruned"] = True
            node["color"] = {
                "background": "rgba(200, 200, 200, 0.5)",
                "border": "rgba(150, 150, 150, 0.5)"
            }
        else:
            node["to_be_pruned"] = False
            node.pop("color", None)

    return Response({
        "network": {
            "nodes": nodes,
            "edges": edges
        },
        "pruned_network": {
            "nodes": [node for node in nodes if node["id"] in pruned_node_ids],
            "edges": pruned_edges
        }
    })


@api_view(["POST"])
def apply_layout(request) -> Response:
    hierachical_layout = request.data.get("hierachical_layout", "False")
    nodes = request.data.get("nodes", "[]")
    if hierachical_layout == "True":
        nodes = generate_hierarchical_layout(nodes)
    else:
        nodes = generate_random_layout(nodes)
    return Response(nodes)


def generate_random_layout(nodes):
    sizing_factor = 30
    G = nx.Graph()
    for node in nodes:
        G.add_node(node["id"])
    pos = nx.random_layout(G, seed=123)
    for node in nodes:
        node["x"] = pos[node["id"]][0] * (len(nodes) * sizing_factor)
        node["y"] = pos[node["id"]][1] * (len(nodes) * sizing_factor)
    return nodes


def generate_hierarchical_layout(nodes):
    sizing_factor = 20
    order_layers = {'Extracellular': 'a', 'Cell surface': 'b', 'Plasma membrane': 'c', 'Cytoplasm': 'd',
                    'Multiple': 'e', 'Nucleus': 'f', 'Other': 'g', 'Unknown': 'h', 'None': 'i'}

    mapper_multiple_layers = {}
    G = nx.Graph()
    for node in nodes:
        if "layer" in node:
            if str(node["layer"]).startswith("Multiple"):
                mapper_multiple_layers[node["id"]] = node["layer"]
                node["layer"] = "Multiple"
            G.add_node(node["id"], layer=order_layers[node["layer"]])
        else:
            G.add_node(node["id"], layer=order_layers["None"])

    pos = nx.multipartite_layout(G, subset_key="layer", align="horizontal", scale=len(nodes) * sizing_factor)

    y_offset = 300
    extra_spacing_layers = {'Unknown', 'Other', 'None'}
    for node in nodes:
        if "id" in node:
            node["x"] = pos[node["id"]][0]
            node["y"] = pos[node["id"]][1]
            if node["id"] in mapper_multiple_layers:
                node["layer"] = mapper_multiple_layers[node["id"]]
            if node.get("layer", "None") in extra_spacing_layers:
                node["y"] += y_offset
    return nodes


@api_view(["POST"])
def map_nodes(request) -> Response:
    """Maps user given input nodes to Proteins in the django database.
    Further updates the node list given by the user by extending the matching proteins with information
    from the database, leaves unmatched nodes untouched. No informations from the input node list gets
    removed. Custom node attributes remain untouched. Returns updated node list.

    Args:
        request (HttpRequest): With keys "nodes" for the node list containing input node objects from the frontend,
        with "id" key, and key "identifier" representing the Protein backend attribute the node id are representing.
        Identifier must be of type "Identifier" as defined in the frontend.

    Returns:
        Response: Updates node list.
    """
    # load data from request
    nodes = request.data.get("nodes", "[]")
    identifier = request.data.get("identifier", "")
    reviewed = request.data.get("reviewed", False)

    nodes = fetch_node_information(nodes, identifier, reviewed)

    # set label to node identifier if label is unset, otherwise
    # return list of nodes updated nodes
    return Response(nodes)


@api_view(["POST"])
def tasks_view(request) -> Response:
    tokens = json.loads(request.data.get("tokens", "[]"))
    tasks = Task.objects.filter(token__in=tokens).order_by("-created_at").all()
    tasks_info = []
    for task in tasks:
        if not task.done and not task.failed:
            refresh_from_redis(task)
            task.save()

        tasks_info.append(
            {
                "token": task.token,
                "info": TaskStatusSerializer().to_representation(task),
                "stats": task_stats(task),
            }
        )
    return Response(tasks_info)


@api_view(["POST"])
def add_edges(request) -> Response:
    if "network" not in request.data:
        return Response(None)
    result = json.loads(request.data["result"])
    parameters = result.get("parameters", {})
    background_mapping = result.get("backgroundMapping", {})
    background_mapping_reverse = result.get("backgroundMappingReverse", {})
    id_space = parameters["config"].get("identifier", "symbol")
    edges = request.data["network"]["edges"]
    nodes = request.data["network"]["nodes"]
    custom_edges = parameters.get("customEdges", False)
    ppi_dataset = parameters.get("ppiDataset")
    pdi_dataset = parameters.get("pdiDataset")
    filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset['licenced'] or pdi_dataset['licenced']:
        filename += "_licenced"
    if parameters["config"].get("reviewed", False):
        filename += "_reviewed"
    path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    data_dir = os.path.join(path, "data", "Networks")
    filename = os.path.join(data_dir, filename + ".gt")
    g = gt.load_graph(filename)
    if custom_edges:
        g = add_edges(g, edges)
    all_nodes_int = set([int(background_mapping[gene["id"]]) for gene in nodes if gene["id"] in background_mapping])
    edges_unique = set()
    for node in nodes:
        for neighbor in g.get_all_neighbors(background_mapping[node["id"]]):
            if int(neighbor) > int(background_mapping[node["id"]]) and int(neighbor) in all_nodes_int:
                first_key = next(iter(background_mapping_reverse))

                if isinstance(first_key, int):
                    neighbor_key = int(neighbor)
                else:
                    neighbor_key = str(int(neighbor))
                edges_unique.add((node["id"], background_mapping_reverse[neighbor_key]))

    edges = [{"from": source, "to": target} for
             source, target in edges_unique]
    return Response(edges)


@api_view(["POST"])
def create_network(request) -> Response:
    if "network" not in request.data:
        return Response(None)
    else:
        if "nodes" not in request.data["network"]:
            request.data["network"]["nodes"] = []
        if "edges" not in request.data["network"]:
            request.data["network"]["edges"] = []
    if "config" not in request.data:
        request.data["config"] = {}
    if "groups" not in request.data:
        request.data["groups"] = {}

    id = uuid.uuid4().hex
    while True:
        try:
            Network.objects.create(
                id=id,
                nodes=request.data["network"]["nodes"],
                edges=request.data["network"]["edges"],
                config=request.data["config"],
                groups=request.data["groups"],
            )
            break
        except IntegrityError:
            id = uuid.uuid4().hex
    return Response(id)


def latest_datasets(ds):
    dataset_dict = {}
    for d in ds:
        name = d.name + "_" + str(d.licenced)
        if name not in dataset_dict:
            dataset_dict[name] = d
            continue
        if dataset_dict[name].version < d.version:
            dataset_dict[name] = d
    return dataset_dict.values()


def get_or_create_network_file(dataset, dataset_type, fmt, reviewed, params):
    from drugstone.management.commands.make_graphs import get_or_create_ppi_network, get_or_create_pdi_network, \
        get_or_create_pdis_network, get_or_create_drdis_network
    # try:
    # reviewed = bool(params.get("reviewed", "True"))
    # except ValueError:
    #     reviewed = True
    print(f"Reviewed proteins only: {reviewed}")
    match dataset_type:
        case "ppi":
            return get_or_create_ppi_network(dataset, params.get("identifier", None), False, reviewed, fmt)
        case "pdi":
            return get_or_create_pdi_network(dataset, params.get("identifier", None), False, reviewed, fmt)
        case "pdis":
            return get_or_create_pdis_network(dataset, params.get("identifier", None), False, reviewed, fmt)
        case "drdis":
            return get_or_create_drdis_network(dataset, False, fmt)
    return "NIY"


@api_view(["GET"])
def download_network(request) -> Response:
    dataset_name = request.query_params.get("dataset")
    dataset_type = request.query_params.get("dataset_type").lower()
    if "_dataset" in dataset_type:
        dataset_type = dataset_type.replace("_dataset", "")
    dataset = None
    match dataset_type:
        case "ppi":
            dataset = get_ppi_ds(dataset_name, False)
        case "pdi":
            dataset = get_pdi_ds(dataset_name, False)
        case "pdis":
            dataset = get_pdis_ds(dataset_name, False)
        case "drdis":
            dataset = get_drdis_ds(dataset_name, False)

    if dataset is None:
        return Response("Dataset not found", status=404)

    format = request.query_params.get("fmt", "gt")
    fmt_list = ["gt", "graphml", "xml", "dot", "gml"]
    if format not in fmt_list:
        return Response(f"Format not supported: {format}! Choose one of: {fmt_list}", status=400)

    reviewed = "false" != request.query_params.get("reviewed", "True").lower()

    file = get_or_create_network_file(dataset, dataset_type, fmt=format, reviewed = reviewed, params=request.query_params)

    if file is not None:
        response = StreamingHttpResponse(FileWrapper(open(file, 'rb'), 512), content_type=mimetypes.guess_type(file)[0])
        _, file_name = os.path.split(file)
        response['Content-Disposition'] = 'attachment; filename=' + smart_str(file_name)
        response['Content-Length'] = os.path.getsize(file)
        return response
    return Response(
        f"A dataset with the given parameters does either not exist or could not be created. Please check your inputs again or try in a few minutes.",
        status=404)


@api_view(["GET"])
def get_datasets(request) -> Response:
    datasets = {}
    datasets["protein-protein"] = PPIDatasetSerializer(many=True).to_representation(
        latest_datasets(PPIDataset.objects.all())
    )
    datasets["protein-drug"] = PDIDatasetSerializer(many=True).to_representation(
        latest_datasets(PDIDataset.objects.all())
    )
    datasets["protein-disorder"] = PDisDatasetSerializer(many=True).to_representation(
        latest_datasets(PDisDataset.objects.all())
    )
    datasets["drug-disorder"] = DrDisDatasetSerializer(many=True).to_representation(
        latest_datasets(DrDiDataset.objects.all())
    )
    return Response(datasets)


@api_view(["GET"])
def load_network(request) -> Response:
    network = NetworkSerializer().to_representation(
        Network.objects.get(id=request.query_params.get("id"))
    )
    result = {
        "network": {
            "nodes": json.loads(network["nodes"].replace("'", '"')),
            "edges": json.loads(network["edges"].replace("'", '"')),
        },
        "config": json.loads(
            network["config"]
            .replace("'", '"')
            .replace("True", "true")
            .replace("False", "false")
        ),
        "groups": json.loads(
            network["groups"]
            .replace("'", '"')
            .replace("True", "true")
            .replace("False", "false")
        ),
    }
    return Response(result)


@api_view(["GET"])
def get_all_scores_pathway_enrichment(request) -> Response:
    token_str = request.query_params["token"]
    task = Task.objects.get(token=token_str)
    result = task_result(task)
    score_preparations = result["score_preparations"]
    seeds = result["parameters"]["seeds"]
    return Response(get_all_node_scores(score_preparations, seeds))


@api_view(["PUT"])
def calculate_result_for_pathway(request) -> Response:
    token_str = request.query_params["token"]
    task = Task.objects.get(token=token_str)
    result = task_result(task)
    geneset = result["mapGenesetsReverse"][request.query_params["geneset"]]
    pathway = request.query_params["pathway"]
    if "geneset" in result and result["geneset"] == geneset and result["pathway"] == pathway:
        # already parsed
        return Response(result)
    path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    data_dir = os.path.join(path, "data", "Networks")

    df_from_json = pd.read_json(result["filteredDf"], orient='records')
    network, isSeed = parse_pathway(geneset, pathway, df_from_json, task.parameters, data_dir,
                                    result["backgroundMapping"], result["backgroundMappingReverse"],
                                    result["mapGenesets"], result["geneSetsDict"], result["score_preparations"])
    result["network"] = network
    result["geneset"] = request.query_params["geneset"]
    result["pathway"] = pathway
    result["node_attributes"] = {}
    result["node_attributes"]["isSeed"] = isSeed
    update_result(result, token_str)
    return Response("worked!")


@api_view(["POST"])
def update_network(request) -> Response:
    token_str = request.data["token"]
    task = Task.objects.get(token=token_str)
    result = task_result(task)
    result["network"] = request.data["network"]
    if "cutoff" in request.data:
        result["cutoff"] = request.data["cutoff"]
    if "prune_orphan_nodes" in request.data:
        result["prune_orphan_nodes"] = request.data["prune_orphan_nodes"]
    update_result(result, token_str)
    return Response("worked!")


@api_view()
def result_view(request) -> Response:
    node_name_attribute = "drugstone_id"

    view = request.query_params.get("view")
    fmt = request.query_params.get("fmt")
    token_str = request.query_params["token"]
    task = Task.objects.get(token=token_str)
    result = task_result(task)
    if result.get("algorithm") == "pathway_enrichment" or result.get("algorithm") == "louvain_clustering" or result.get(
            "algorithm") == "leiden_clustering" or result.get("algorithm") == "first_neighbor":
        return Response(result)
    node_attributes = result.get("node_attributes")
    if not node_attributes:
        node_attributes = {}
        result["node_attributes"] = node_attributes

    proteins = []
    drugs = []

    network = result["network"]
    node_types = node_attributes.get("node_types")
    if not node_types:
        node_types = {}
        node_attributes["node_types"] = node_types

    is_seed = node_attributes.get("is_seed")
    if not is_seed:
        is_seed = {}
        node_attributes["is_seed"] = is_seed
    scores = node_attributes.get("scores", {})
    node_details = {}
    protein_id_map = defaultdict(set)
    node_attributes["details"] = node_details
    parameters = json.loads(task.parameters)
    seeds = parameters["seeds"]
    nodes = network["nodes"]

    parameters = task_parameters(task)
    # attach input parameters to output
    result["parameters"] = parameters
    identifier_nodes = set()
    identifier = parameters["config"]["identifier"]
    if not "reviewed" in parameters["config"]:
        parameters["config"]["reviewed"] = False

    # merge input network with result network
    for node in parameters["input_network"]["nodes"]:
        # if node was already mapped, add user defined values to result of analysis
        if identifier in node:
            node_name = node[identifier][0]
            if node_name in node_details:
                # update the node to not lose user input attributes
                node_details[node_name].update(node)
                # skip adding node if node already exists in analysis output to avoid duplicates
            else:
                # node does not exist in analysis output yet, was added by user but not used as seed
                node_details[node_name] = node
                # append mapped input node to analysis result
                nodes.append(node_name)
                # manually add node to node types
                result["node_attributes"]["node_types"][node_name] = "protein"
        else:
            # node is custom node from user, not mapped to drugstone but will be displayed with all custom attributes
            node_id = node["id"]
            identifier_nodes.add(node_id)
            node_details[node_id] = node
            is_seed[node_id] = False
            # append custom node to analysis result later on
            # manually add node to node types
            result["node_attributes"]["node_types"][node_id] = "custom"
    # extend the analysis network by the input netword nodes
    # map edge endpoints to database proteins if possible and add edges to analysis network
    protein_nodes = set()
    # mapping all new protein and drug nodes by drugstoneIDs + adding scores
    for node_id in nodes:
        if node_id[:2] == "dr":
            node_data = DrugSerializer().to_representation(
                Drug.objects.get(id=int(node_id[2:]))
            )
            node_data["drugstoneType"] = "drug"
            drugs.append(node_data)
            if node_id in scores:
                node_data["score"] = scores.get(node_id, None)
            node_types[node_id] = "drug"
            node_details[node_id] = node_data
        elif node_id[:2] != "di":
            protein_nodes.add(node_id)
        else:
            continue

    nodes_mapped, identifier = query_proteins_by_identifier(protein_nodes, identifier, parameters["config"]["reviewed"])

    nodes_mapped_dict = {node[identifier][0]: node for node in nodes_mapped}

    # merge fetched data with given data to avoid data loss
    for node_id in nodes:
        if node_id in nodes_mapped_dict:
            # node.update(nodes_mapped_dict[node['id']])
            node_data = nodes_mapped_dict[node_id]
            node_data["drugstoneType"] = "protein"
            # proteins.append(node_data)
            node_ident = node_data[identifier][0]
            # node_data[identifier] = [node_ident]
            protein_id_map[node_ident].add(node_id)
            identifier_nodes.add(node_ident)
            is_seed[node_ident] = node_id in seeds or (
                is_seed[node_ident] if node_ident in is_seed else False
            )
            node_types[node_ident] = "protein"
            score = scores.get(node_id, None)
            if node_ident in node_details:
                data = node_details[node_ident]
                data["score"] = [score] if score else None
            else:
                node_data["score"] = score if score else None
                node_data["drugstoneType"] = "protein"
                node_data["id"] = node_ident
                node_data["label"] = node_ident
                node_details[node_ident] = node_data

    for node_id, detail in node_details.items():
        if "drugstoneType" in detail and detail["drugstoneType"] == "protein":
            detail["symbol"] = list(set(detail["symbol"])) if "symbol" in detail else []
            detail["entrez"] = list(set(detail["entrez"])) if "entrez" in detail else []
            detail["uniprot"] = (
                list(set(detail["uniprot"])) if "uniprot" in detail else []
            )
            detail["ensg"] = list(set(detail["ensg"])) if "ensg" in detail else []

    edges = parameters["input_network"]["edges"]

    edge_endpoint_ids = set()

    # TODO check for custom edges when working again with ensemble gene ids
    for edge in edges:
        edge_endpoint_ids.add(edge["from"])
        edge_endpoint_ids.add(edge["to"])

    nodes_mapped, id_key = query_proteins_by_identifier(edge_endpoint_ids, identifier, parameters["config"]["reviewed"])

    pdi_config = result.get("parameters").get('pdi_dataset')

    if pdi_config:
        pdi_dataset = get_pdi_ds(pdi_config.get('name', DEFAULTS['pdi']), pdi_config.get('licenced', False))
        for edge in result['network']['edges']:
            if (edge['from'][:2] == 'dr'):
                # drug should always be "to", flip edge
                drug = edge['from']
                edge['from'] = edge['to']
                edge['to'] = drug
            if (edge['to'][:2] == 'dr'):
                drug_id = int(edge['to'][2:])
                pdi_object = ProteinDrugInteraction.objects.filter(
                    protein_id__in={int(p[1:]) for p in node_attributes['details'][edge['from']]['drugstone_id']},
                    drug_id=drug_id, pdi_dataset_id=pdi_dataset.id)
                actions = set()
                for pdi in pdi_object:
                    if pdi.actions:
                        for action in json.loads(pdi.actions):
                            actions.add(action)
                edge['actions'] = list(actions)

    if (
            "autofill_edges" in parameters["config"]
            and parameters["config"]["autofill_edges"]
    ):
        prots = list(
            filter(
                lambda n: n["drugstone_type"] == "protein",
                filter(
                    lambda n: "drugstone_type" in n and node_name_attribute in n,
                    parameters["input_network"]["nodes"],
                ),
            )
        )
        proteins = {
            node_name[1:] for node in prots for node_name in node[node_name_attribute]
        }
        dataset = (
            DEFAULTS["ppi"]
            if "interaction_protein_protein" not in parameters["config"]
            else parameters["config"]["interaction_protein_protein"]
        )
        dataset_object = models.PPIDataset.objects.filter(name__iexact=dataset).last()
        interaction_objects = models.ProteinProteinInteraction.objects.filter(
            Q(ppi_dataset=dataset_object)
            & Q(from_protein__in=proteins)
            & Q(to_protein__in=proteins)
        )
        auto_edges = list(
            map(
                lambda n: {
                    "from": f"p{n.from_protein_id}",
                    "to": f"p{n.to_protein_id}",
                    "is_directed": f"{n.is_directed}",
                    "is_stimulation": f"{n.is_stimulation}",
                    "is_inhibition": f"{n.is_inhibition}",
                },
                interaction_objects,
            )
        )
        edges.extend(auto_edges)

    result["network"]["edges"].extend(edges)
    uniq_edges = dict()
    for edge in result["network"]["edges"]:
        hash = edge["from"] + edge["to"]
        uniq_edges[hash] = edge
    result["network"]["edges"] = list(uniq_edges.values())

    # Only map edges if the edge source is Omnipath (directed)
    if result.get("parameters").get('ppi_dataset')['name'] == "OmniPath":
        drugstone_edges = []
        for edge in result["network"]["edges"]:
            if edge["from"] in nodes_mapped_dict and edge["to"] in nodes_mapped_dict:
                fr = nodes_mapped_dict[edge["from"]]['drugstone_id'][0]
                to = nodes_mapped_dict[edge["to"]]['drugstone_id'][0]
                edge_data = {k: v for k, v in edge.items() if k not in ['from', 'to']}
                edge_data.update({"from": fr, "to": to})
                drugstone_edges.append(edge_data)
            else:
                drugstone_edges.append(edge)
        result["network"]["edges"] = fetch_edges_from_input(result.get("parameters").get('ppi_dataset')['name'],
                                                            result.get("parameters").get('ppi_dataset')['licenced'],
                                                            drugstone_edges)

    if "scores" in result["node_attributes"]:
        del result["node_attributes"]["scores"]

    if "properties" in result:
        for node in result["node_attributes"]["details"].values():
            if "id" in node.keys() and node["id"] in result["properties"]:
                if "properties" not in node:
                    node["properties"] = {}
                node["properties"].update(result["properties"][node["id"]])
    if not view:
        return Response(result)
    else:
        if view == "proteins":
            proteins = list(
                filter(
                    lambda n: "drugstone_type" in n
                              and n["drugstone_type"] == "protein",
                    node_details.values(),
                )
            )
            if fmt == "csv":
                items = []
                for i in proteins:
                    new_i = {
                        "id": i["id"],
                        "uniprot": i["uniprot"] if "uniprot" in i else [],
                        "gene": i["symbol"] if "symbol" in i else [],
                        "name": i["protein_name"] if "protein_name" in i else [],
                        "ensembl": i["ensg"] if "ensg" in i else [],
                        "entrez": i["entrez"] if "entrez" in i else [],
                        "seed": is_seed[i["id"]],
                    }
                    if "score" in i:
                        new_i["score"] = i["score"]
                    items.append(new_i)
            else:
                items = proteins
        elif view == "drugs":
            if fmt == "csv":
                items = [i for i in drugs]
            else:
                items = drugs
        else:
            return Response({})

        if not fmt or fmt == "json":
            return Response(items)
        elif fmt == "csv":
            if len(items) != 0:
                keys = items[0].keys()
            else:
                keys = []
            response = HttpResponse(content_type="text/csv")
            response[
                "Content-Disposition"
            ] = f'attachment; filename="{task.token}_{view}.csv"'
            dict_writer = csv.DictWriter(response, keys)
            dict_writer.writeheader()
            dict_writer.writerows(items)
            return response
        else:
            return Response({})


@api_view(["POST"])
def graph_export(request) -> Response:
    """
    Recieve whole graph data and write it to graphml file. Return the
    file ready to download.
    """
    remove_node_properties = [
        "color",
        "shape",
        "border_width",
        "group",
        "border_width_selected",
        "shadow",
        "group_id",
        "drugstone_type",
        "font",
        "x",
        "y",
        "_group",
    ]
    rename_node_properties = {"group_name": "group"}
    remove_edge_properties = ["group", "color", "dashes", "shadow", "id"]
    rename_edge_properties = {"group_name": "group"}
    nodes = request.data.get("nodes", [])
    edges = request.data.get("edges", [])
    fmt = request.data.get("fmt", "graphml")
    G = nx.Graph()
    node_map = dict()
    for node in nodes:
        # networkx does not support datatypes such as lists or dicts
        for prop in remove_node_properties:
            if prop in node:
                del node[prop]
        for k, v in rename_node_properties.items():
            if k in node:
                node[v] = node[k]
                del node[k]
        for key in list(node.keys()):
            if isinstance(node[key], list) or isinstance(node[key], dict):
                node[key] = json.dumps(node[key])
            elif node[key] is None:
                # networkx has difficulties with None when writing graphml
                node[key] = ""
        try:
            node_name = node["label"]
            if "drugstone_id" in node:
                node_map[node["drugstone_id"]] = node["label"]
            elif "id" in node:
                node_map[node["id"]] = node["label"]
        except KeyError:
            node_name = node["drugstone_id"]
        G.add_node(node_name, **node)

    for e in edges:
        # networkx does not support datatypes such as lists or dicts
        for prop in remove_edge_properties:
            if prop in e:
                del e[prop]
        for k, v in rename_edge_properties.items():
            if k in e:
                e[v] = e[k]
                del e[k]
        for key in e:
            if isinstance(e[key], list) or isinstance(e[key], dict):
                e[key] = json.dumps(e[key])
            elif e[key] is None:
                e[key] = ""
        u_of_edge = e.pop("from")
        u_of_edge = u_of_edge if u_of_edge not in node_map else node_map[u_of_edge]
        v_of_edge = e.pop("to")
        v_of_edge = node_map[v_of_edge] if v_of_edge in node_map else v_of_edge
        G.add_edge(u_of_edge, v_of_edge, **e)

    if fmt == "graphml":
        data = nx.generate_graphml(G)
        response = HttpResponse(data, content_type="application/xml")
    elif fmt == "json":
        data = nx.readwrite.json_graph.node_link_data(G)
        del data["graph"]
        del data["multigraph"]

        # for node in data['nodes']:
        # for prop in remove_node_properties:
        #     if prop in node:
        #         del node[prop]
        # for edge in data['links']:
        # for prop in remove_edge_properties:
        #     if prop in edge:
        #         del edge[prop]
        data["edges"] = data.pop("links")
        data = json.dumps(data)
        data = (
            data.replace('"{', "{")
            .replace('}"', "}")
            .replace('"[', "[")
            .replace(']"', "]")
            .replace('\\"', '"')
        )
        response = HttpResponse(data, content_type="application/json")
    elif fmt == "csv":
        data = pd.DataFrame(
            nx.to_numpy_array(G), columns=G.nodes(), index=G.nodes(), dtype=int
        )
        response = HttpResponse(data.to_csv(), content_type="text/csv")

    response[
        "content-disposition"
    ] = f'attachment; filename="{int(time.time())}_network.{fmt}"'
    return response


@api_view(["POST"])
def adjacent_disorders(request) -> Response:
    """Find all adjacent disorders to a list of proteins.

    Args:
        request (django.request): Request object with keys "proteins" and "pdi_dataset"

    Returns:
        Response: With lists "pdis" (protein-drug-intersions) and "disorders"
    """
    data = request.data
    if "proteins" in data:
        drugstone_ids = data.get("proteins", [])
        pdis_dataset = get_pdis_ds(
            data.get("dataset", DEFAULTS["pdis"]), data.get("licenced", False)
        )
        # find adjacent drugs by looking at drug-protein edges
        pdis_objects = ProteinDisorderAssociation.objects.filter(
            protein__id__in=drugstone_ids, pdis_dataset_id=pdis_dataset.id
        )
        disorders = {e.disorder for e in pdis_objects}
        # serialize
        edges = ProteinDisorderAssociationSerializer(many=True).to_representation(
            pdis_objects
        )
        disorders = DisorderSerializer(many=True).to_representation(disorders)
    elif "drugs" in data:
        drugstone_ids = data.get("drugs", [])
        drdi_dataset = get_drdis_ds(
            data.get("dataset", DEFAULTS["drdi"]), data.get("licenced", False)
        )
        # find adjacent drugs by looking at drug-protein edges
        drdi_objects = DrugDisorderIndication.objects.filter(
            drug__id__in=drugstone_ids, drdi_dataset_id=drdi_dataset.id
        )
        disorders = {e.disorder for e in drdi_objects}
        # serialize
        edges = DrugDisorderIndicationSerializer(many=True).to_representation(
            drdi_objects
        )
        disorders = DisorderSerializer(many=True).to_representation(disorders)
    for d in disorders:
        d["drugstone_type"] = "disorder"
    return Response(
        {
            "edges": edges,
            "disorders": disorders,
        }
    )


@api_view(["POST"])
def adjacent_drugs(request) -> Response:
    """Find all adjacent drugs to a list of proteins.

    Args:
        request (django.request): Request object with keys "proteins" and "pdi_dataset"

    Returns:
        Response: With lists "pdis" (protein-drug-intersions) and "drugs"
    """
    data = request.data
    drugstone_ids = data.get("proteins", [])
    pdi_dataset = get_pdi_ds(
        data.get("pdi_dataset", DEFAULTS["pdi"]), data.get("licenced", False)
    )
    approved = data.get("approved", False)
    # find adjacent drugs by looking at drug-protein edges
    pdi_objects = ProteinDrugInteraction.objects.filter(
        protein__id__in=drugstone_ids, pdi_dataset_id=pdi_dataset.id
    )
    drugs = {e.drug for e in pdi_objects}
    # serialize
    pdis = ProteinDrugInteractionSerializer(many=True).to_representation(pdi_objects)
    drugs = DrugSerializer(many=True).to_representation(drugs)
    if approved:
        drugs = [drug for drug in drugs if drug["status"] == "approved"]
    for drug in drugs:
        drug["drugstone_type"] = "drug"

    return Response(
        {
            "pdis": pdis,
            "drugs": drugs,
        }
    )


@api_view(["POST"])
def query_proteins(request) -> Response:
    proteins = request.data

    details = []
    not_found = []
    for p in proteins:
        try:
            protein = Protein.objects.get(uniprot_code=p)
            details.append(ProteinSerializer().to_representation(protein))
            continue
        except Protein.DoesNotExist:
            pass

        drug_interactions = ProteinDrugInteraction.objects.filter(drug__drug_id=p)
        if len(drug_interactions) > 0:
            for di in drug_interactions:
                details.append(ProteinSerializer().to_representation(di.protein))
            continue

        not_found.append(p)

    return Response(
        {
            "details": details,
            "notFound": not_found,
        }
    )


@api_view(["POST"])
def send_bugreport(request) -> Response:
    data = request.data
    title = data.get("title")
    body = data.get("body")
    email = data.get("email", None)
    if email and len(email) == 0:
        email = None
    if not title or not body:
        return Response({"status": 400})

    bugreport(title, body, email)
    return Response({"status": 200})


@api_view(["POST"])
def save_selection(request) -> Response:
    chars = string.ascii_lowercase + string.ascii_uppercase + string.digits
    token_str = "".join(random.choice(chars) for _ in range(32))

    config = request.data.get("config")
    network = request.data.get("network")

    Network.objects.create(id=token_str, config=json.dumps(config), nodes=json.dumps(network["nodes"]),
                           edges=json.dumps(network["edges"]))
    return Response({
        'token': token_str,
    })


@api_view(["PUT"])
def rename_selection(request) -> Response:
    print(request.data)
    token = request.data.get("token")
    name = request.data.get("name")

    if not token or not name:
        return Response({"error": "Missing 'token' or 'name'"}, status=400)

    try:
        network = Network.objects.get(id=token)
        network.name = name
        network.save()
        return Response({"message": "Name updated successfully."})
    except Network.DoesNotExist:
        return Response({"error": "Network not found"}, status=404)


@api_view(["GET"])
def get_view(request) -> Response:
    token = request.query_params.get("token")
    network = Network.objects.get(id=token)
    return Response(
        {
            "config": json.loads(network.config),
            "created_at": network.created_at,
            "name": network.name,
            "network": {
                "nodes": json.loads(network.nodes),
                "edges": json.loads(network.edges),
            },
        }
    )


@api_view(["POST"])
def get_view_infos(request) -> Response:
    tokens = request.data.get('tokens')
    networks = Network.objects.filter(id__in=tokens).order_by('-created_at')
    return Response([{
        'token': n.id,
        'created_at': n.created_at,
        'name': n.name,
    } for n in networks])


@api_view(["GET"])
def get_max_tissue_expression(request) -> Response:
    tissue = Tissue.objects.get(id=request.query_params.get("tissue"))
    return Response(
        {
            "max": ExpressionLevel.objects.filter(tissue=tissue).aggregate(
                Max("expression_level")
            )["expression_level__max"]
        }
    )


@api_view(["POST"])
def query_tissue_proteins(request) -> Response:
    threshold = request.data["threshold"]
    tissue_id = request.data["tissue_id"]
    tissue = Tissue.objects.get(id=tissue_id)

    proteins = []
    for el in tissue.expressionlevel_set.filter(expression_level__gte=threshold):
        proteins.append(ProteinSerializer().to_representation(el.protein))

    return Response(proteins)


class TissueView(APIView):
    def get(self, request) -> Response:
        tissues = Tissue.objects.all()
        return Response(TissueSerializer(many=True).to_representation(tissues))


class TissueExpressionView(APIView):
    """
    Expression of host proteins in tissues.
    """

    def get(self, request) -> Response:
        tissue = Tissue.objects.get(id=request.query_params.get("tissue"))
        proteins = request.query_params.get("proteins")
        token = request.query_params.get("token")
        return self.get_tissue_expression(tissue, proteins, token)

    def post(self, request) -> Response:
        tissue = Tissue.objects.get(id=request.data.get("tissue"))
        proteins = request.data.get("proteins")
        token = request.data.get("token")
        return self.get_tissue_expression(tissue, proteins, token)

    def get_tissue_expression(self, tissue, proteins, token):
        if proteins is not None:
            ids = json.loads(proteins)
            proteins = list(Protein.objects.filter(id__in=ids).all())
        elif token is not None:
            proteins = []
            task = Task.objects.get(token=token)
            result = task_result(task)
            network = result["network"]
            node_attributes = result.get("node_attributes")
            if not node_attributes:
                node_attributes = {}
            node_types = node_attributes.get("node_types")
            if not node_types:
                node_types = {}
            parameters = json.loads(task.parameters)
            seeds = parameters["seeds"]
            nodes = network["nodes"]
            for node in nodes + seeds:
                node_type = node_types.get(node)
                details = None
                if node_type == "protein":
                    if details:
                        proteins.append(details)
                    else:
                        try:
                            prot = Protein.objects.get(uniprot_code=node)
                            if prot not in proteins:
                                proteins.append(Protein.objects.get(uniprot_code=node))
                        except Protein.DoesNotExist:
                            pass

        pt_expressions = {}

        for protein in proteins:
            try:
                expression_level = ExpressionLevel.objects.get(
                    protein=protein, tissue=tissue
                )
                pt_expressions[
                    ProteinSerializer().to_representation(protein)["drugstone_id"]
                ] = expression_level.expression_level
            except ExpressionLevel.DoesNotExist:
                pt_expressions[
                    ProteinSerializer().to_representation(protein)["drugstone_id"]
                ] = None

        return Response(pt_expressions)
