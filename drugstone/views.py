import csv
import random
import string
import time
import uuid
from collections import defaultdict
import pandas as pd
import networkx as nx
from django.http import HttpResponse
from django.db.models import Q, Max
from django.db import IntegrityError
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.views import APIView

from drugstone.util.mailer import bugreport
from drugstone.util.query_db import (
    query_proteins_by_identifier,
    clean_proteins_from_compact_notation,
    fetch_node_information
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

from drugstone.settings import DEFAULTS


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


@api_view(["POST"])
def convert_compact_ids(request) -> Response:
    nodes = request.data.get("nodes", "[]")
    identifier = request.data.get("identifier", "")
    cleaned = clean_proteins_from_compact_notation(nodes, identifier)
    return Response(cleaned)


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
    
    nodes = fetch_node_information(nodes, identifier)

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


@api_view()
def result_view(request) -> Response:
    node_name_attribute = "drugstone_id"

    view = request.query_params.get("view")
    fmt = request.query_params.get("fmt")
    token_str = request.query_params["token"]
    task = Task.objects.get(token=token_str)
    result = task_result(task)
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

    nodes_mapped, identifier = query_proteins_by_identifier(protein_nodes, identifier)

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

    nodes_mapped, id_key = query_proteins_by_identifier(edge_endpoint_ids, identifier)

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

    if "scores" in result["node_attributes"]:
        del result["node_attributes"]["scores"]

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
    # find adjacent drugs by looking at drug-protein edges
    pdi_objects = ProteinDrugInteraction.objects.filter(
        protein__id__in=drugstone_ids, pdi_dataset_id=pdi_dataset.id
    )
    drugs = {e.drug for e in pdi_objects}
    # serialize
    pdis = ProteinDrugInteractionSerializer(many=True).to_representation(pdi_objects)
    drugs = DrugSerializer(many=True).to_representation(drugs)
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


@api_view(["GET"])
def get_view(request) -> Response:
    token = request.query_params.get("token")
    network = Network.objects.get(id=token)
    return Response(
        {
            "config": json.loads(network.config),
            "created_at": network.created_at,
            "network": {
                "nodes": json.loads(network.nodes),
                "edges": json.loads(network.edges),
            },
        }
    )


@api_view(["POST"])
def get_view_infos(request) -> Response:
    tokens = request.data.get('tokens')
    networks = Network.objects.filter(id__in=tokens)
    return Response([{
        'token': n.id,
        'created_at': n.created_at,
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
