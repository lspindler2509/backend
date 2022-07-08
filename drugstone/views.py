import csv
import json
import random
import string
import time
import uuid
import pandas as pd
from typing import Tuple

import networkx as nx
from django.http import HttpResponse
from django.db.models import Q
from django.db import IntegrityError
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework.views import APIView
from drugstone.util.query_db import query_proteins_by_identifier

from drugstone import models
from drugstone import serializers

from drugstone.models import Protein, Task, ProteinDrugInteraction, \
    Drug, Tissue, ExpressionLevel, Network, ProteinDisorderAssociation, DrugDisorderIndication, Disorder, DrDiDataset, PDIDataset, PDisDataset, PPIDataset
from drugstone.serializers import ProteinSerializer, TaskSerializer, \
    ProteinDrugInteractionSerializer, DrugSerializer, TaskStatusSerializer, TissueSerializer, NetworkSerializer, \
    ProteinDisorderAssociationSerializer, DisorderSerializer, DrugDisorderIndicationSerializer
from drugstone.backend_tasks import start_task, refresh_from_redis, task_stats, task_result, task_parameters


# we might want to replace this class with some ProteinProteinInteraction view of user input proteins

# class ProteinViralInteractionView(APIView):
#     """
#     Protein-Virus-Interaction Network
#     """
#
#     def get(self, request):
#         if not request.query_params.get('data'):
#             proteins = Protein.objects.all()
#             effects = ViralProtein.objects.all()
#             edges = ProteinViralInteraction.objects.all()
#
#             network = {
#                 'proteins': ProteinSerializer(many=True).to_representation(proteins),
#                 'effects': ViralProteinSerializer(many=True).to_representation(effects),
#                 'edges': ProteinViralInteractionSerializer(many=True).to_representation(edges),
#             }
#             return Response(network)
#
#         dataset_virus_list = json.loads(request.query_params.get('data', '[]'))
#         effects = []
#         for dataset_name, virus_name in dataset_virus_list:
#             dataset_virus_object = DatasetVirus.objects.get(dataset=dataset_name, virus=virus_name)
#             effects.extend(list(ViralProtein.objects.filter(dataset_virus=dataset_virus_object).all()))
#
#         edges = []
#         proteins = []
#         for effect in effects:
#             edge_objects = ProteinViralInteraction.objects.filter(effect=effect)
#             for edge_object in edge_objects:
#                 edges.append(edge_object)
#
#                 if edge_object.protein not in proteins:
#                     proteins.append(edge_object.protein)
#
#         network = {
#             'proteins': ProteinSerializer(many=True).to_representation(proteins),
#             'effects': ViralProteinSerializer(many=True).to_representation(effects),
#             'edges': ProteinViralInteractionSerializer(many=True).to_representation(edges),
#         }
#         return Response(network)


# class ProteinDrugInteractionView(APIView):
#     """
#     Protein-Drug-Interaction Network
#     """
#
#     def get(self, request) -> Response:
#         if request.query_params.get('proteins'):
#             print("getting drugs for proteins")
#             protein_ac_list = json.loads(request.query_params.get('proteins'))
#             proteins = list(Protein.objects.filter(uniprot_code__in=protein_ac_list).all())
#         else:
#             proteins = []
#             task = Task.objects.get(token=request.query_params['token'])
#             result = task_result(task)
#             network = result['network']
#             node_attributes = result.get('node_attributes')
#             if not node_attributes:
#                 node_attributes = {}
#             node_types = node_attributes.get('node_types')
#             if not node_types:
#                 node_types = {}
#             nodes = network['nodes']
#             for node in nodes:
#                 node_type = node_types.get(node)
#                 details = None
#                 # if not node_type:
#                 #     print('we should not see this 1')
#                 #     node_type, details = infer_node_type_and_details(node)
#                 if node_type == 'protein':
#                     if details:
#                         proteins.append(details)
#                     else:
#                         try:
#                             proteins.append(Protein.objects.get(uniprot_code=node))
#                         except Protein.DoesNotExist:
#                             pass
#
#         pd_interactions = []
#         drugs = []
#
#         for protein in proteins:
#             pdi_object_list = ProteinDrugInteraction.objects.filter(protein=protein)
#             for pdi_object in pdi_object_list:
#                 pd_interactions.append(pdi_object)
#                 drug = pdi_object.drug
#                 if drug not in drugs:
#                     drugs.append(drug)
#
#         protein_drug_edges = {
#             'proteins': ProteinSerializer(many=True).to_representation(proteins),
#             'drugs': DrugSerializer(many=True).to_representation(drugs),
#             'edges': ProteinDrugInteractionSerializer(many=True).to_representation(pd_interactions),
#         }
#         return Response(protein_drug_edges)


class TaskView(APIView):

    def post(self, request) -> Response:
        chars = string.ascii_lowercase + string.ascii_uppercase + string.digits
        token_str = ''.join(random.choice(chars) for _ in range(32))
        parameters = request.data['parameters']

        # find databases based on parameter strings
        parameters['ppi_dataset'] = serializers.PPIDatasetSerializer().to_representation(
            models.PPIDataset.objects.filter(name__iexact=parameters.get('ppi_dataset', 'STRING')).last())
        parameters['pdi_dataset'] = serializers.PDIDatasetSerializer().to_representation(
            models.PDIDataset.objects.filter(name__iexact=parameters.get('pdi_dataset', 'DrugBank')).last())

        task = Task.objects.create(token=token_str,
                                   target=request.data['target'],
                                   algorithm=request.data['algorithm'],
                                   parameters=json.dumps(parameters))
        start_task(task)
        task.save()

        return Response({
            'token': token_str,
        })

    def get(self, request) -> Response:
        token_str = request.query_params['token']
        task = Task.objects.get(token=token_str)

        if not task.done and not task.failed:
            refresh_from_redis(task)
            task.save()

        return Response({
            'token': task.token,
            'info': TaskSerializer().to_representation(task),
            'stats': task_stats(task),
        })


@api_view(['POST'])
def fetch_edges(request) -> Response:
    """Retrieves interactions between nodes given as a list of drugstone IDs.

    Args:
        request (HttpRequest): With keys 'nodes' representing nodes and 'dataset' representing the
        protein-protein interaction dataset.

    Returns:
        Response: List of edges which are objects with 'from' and to ' attribtues'
    """
    dataset = request.data.get('dataset', 'STRING')
    drugstone_ids = [node['drugstone_id'][1:] for node in request.data.get('nodes', '[]') if 'drugstone_id' in node]
    dataset_object = models.PPIDataset.objects.filter(name__iexact=dataset).last()
    interaction_objects = models.ProteinProteinInteraction.objects.filter(
        Q(ppi_dataset=dataset_object) & Q(from_protein__in=drugstone_ids) & Q(to_protein__in=drugstone_ids))

    return Response(serializers.ProteinProteinInteractionSerializer(many=True).to_representation(interaction_objects))


@api_view(['POST'])
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
    nodes = request.data.get('nodes', '[]')
    id_map = {}
    for node in nodes:
        upper = node['id'].upper()
        id_map[upper] = node['id']
        node['id'] = upper

    identifier = request.data.get('identifier', '')
    # extract ids for filtering
    node_ids = set([node['id'] for node in nodes])

    # query protein table
    nodes_mapped, id_key = query_proteins_by_identifier(node_ids, identifier)

    # change data structure to dict in order to be quicker when merging
    if identifier == 'ensg':
        # a protein might have multiple ensg-numbers, unpack these into single nodes
        nodes_mapped_dict = {node_id: node for node in nodes_mapped for node_id in node[id_key]}
    else:
        nodes_mapped_dict = {node[id_key]: node for node in nodes_mapped}
    # merge fetched data with given data to avoid data loss
    for node in nodes:
        if node['id'] in nodes_mapped_dict:
            node.update(nodes_mapped_dict[node['id']])
        node['id'] = id_map[node['id']]
    # set label to node identifier if label is unset, otherwise
    # return list of nodes updated nodes
    return Response(nodes)


@api_view(['POST'])
def tasks_view(request) -> Response:
    tokens = json.loads(request.data.get('tokens', '[]'))
    tasks = Task.objects.filter(token__in=tokens).order_by('-created_at').all()
    tasks_info = []
    for task in tasks:
        if not task.done and not task.failed:
            refresh_from_redis(task)
            task.save()

        tasks_info.append({
            'token': task.token,
            'info': TaskStatusSerializer().to_representation(task),
            'stats': task_stats(task),
        })
    return Response(tasks_info)


# def infer_node_type_and_details(node) -> Tuple[str, Protein or Drug]:
#     node_type_indicator = node[0]
#     if node_type_indicator == 'p':
#         node_id = int(node[1:])
#         # protein
#         prot = Protein.objects.get(id=node_id)
#         return 'protein', prot
#     elif node_type_indicator == 'd':
#         node_id = int(node[2:])
#         # drug
#         if node_id[0] == 'r':
#             drug = Drug.objects.get(id=node_id[1:])
#             return 'drug', drug
#         elif node_id[0] == 'i':
#             disorder = Disorder.objects.get(id=node_id[1:])
#             return 'disorder', disorder
#     return None, None


@api_view(['POST'])
def create_network(request) -> Response:
    if 'network' not in request.data:
        return Response(None)
    else:
        if 'nodes' not in request.data['network']:
            request.data['network']['nodes'] = []
        if 'edges' not in request.data['network']:
            request.data['network']['edges'] = []
    if 'config' not in request.data:
        request.data['config'] = {}

    id = uuid.uuid4().hex
    while True:
        try:
            Network.objects.create(id=id, nodes=request.data['network']['nodes'],
                                   edges=request.data['network']['edges'], config=request.data['config'])
            break
        except IntegrityError:
            id = uuid.uuid4().hex
    return Response(id)


@api_view(['GET'])
def load_network(request) -> Response:
    network = NetworkSerializer().to_representation(Network.objects.get(id=request.query_params.get('id')))
    result = {'network': {'nodes': json.loads(network['nodes'].replace("'", '"')),
                          'edges': json.loads(network['edges'].replace("'", '"'))},
              'config': json.loads(
                  network['config'].replace("'", '"').replace('True', 'true').replace('False', 'false'))}
    return Response(result)


@api_view()
def result_view(request) -> Response:
    node_name_attribute = 'drugstone_id'

    view = request.query_params.get('view')
    fmt = request.query_params.get('fmt')
    token_str = request.query_params['token']
    task = Task.objects.get(token=token_str)
    result = task_result(task)
    node_attributes = result.get('node_attributes')
    if not node_attributes:
        node_attributes = {}
        result['node_attributes'] = node_attributes
    proteins = []
    drugs = []

    network = result['network']
    node_types = node_attributes.get('node_types')
    if not node_types:
        node_types = {}
        node_attributes['node_types'] = node_types
    is_seed = node_attributes.get('is_seed')
    if not is_seed:
        is_seed = {}
        node_attributes['is_seed'] = is_seed
    scores = node_attributes.get('scores', {})
    node_details = {}
    node_attributes['details'] = node_details
    parameters = json.loads(task.parameters)
    seeds = parameters['seeds']
    nodes = network['nodes']
    # edges = network['edges']
    for node_id in nodes:
        is_seed[node_id] = node_id in seeds
        node_type = node_types.get(node_id).lower()
        pvd_entity = None
        details_s = None
        if node_type == 'protein':
            pvd_entity = Protein.objects.get(id=int(node_id[1:]))
        elif node_type == 'drug':
            pvd_entity = Drug.objects.get(id=int(node_id[2:]))

        if not node_type or not pvd_entity:
            continue
        if node_type == 'protein':
            details_s = ProteinSerializer().to_representation(pvd_entity)
        elif node_type == 'drug':
            details_s = DrugSerializer().to_representation(pvd_entity)
        node_types[node_id] = node_type

        if scores.get(node_id) is not None:
            details_s['score'] = scores.get(node_id, None)
        node_details[node_id] = details_s
        if node_type == 'protein':
            proteins.append(details_s)
        elif node_type == 'drug':
            drugs.append(details_s)
    parameters = task_parameters(task)
    # attach input parameters to output
    result['parameters'] = parameters

    # TODO move the merging to "scores to result"
    # merge input network with result network
    for node in parameters['input_network']['nodes']:
        # if node was already mapped, add user defined values to result of analysis
        if node_name_attribute in node:
            if node[node_name_attribute] in node_details:
                # update the node to not lose user input attributes
                node_details[node[node_name_attribute]].update(node)
                # skip adding node if node already exists in analysis output to avoid duplicates
            else:
                # node does not exist in analysis output yet, was added by user but not used as seed
                node_details[node[node_name_attribute]] = node
                # append mapped input node to analysis result
                nodes.append(node[node_name_attribute])
                # manually add node to node types
                result['node_attributes']['node_types'][node[node_name_attribute]] = 'protein'
        else:
            # node is custom node from user, not mapped to drugstone but will be displayed with all custom attributes
            node_id = node['id']
            nodes.append(node_id)
            node_details[node_id] = node
            is_seed[node_id] = False
            # append custom node to analysis result later on
            # manually add node to node types
            result['node_attributes']['node_types'][node_id] = 'custom'
    # extend the analysis network by the input netword nodes
    # map edge endpoints to database proteins if possible and add edges to analysis network
    identifier = parameters['config']['identifier']
    edges = parameters['input_network']['edges']
    edge_endpoint_ids = set()
    for edge in edges:
        edge_endpoint_ids.add(edge['from'])
        edge_endpoint_ids.add(edge['to'])

    # query protein table
    nodes_mapped, id_key = query_proteins_by_identifier(edge_endpoint_ids, identifier)
    # change data structure to dict in order to be quicker when merging
    nodes_mapped_dict = {node[id_key]: node for node in nodes_mapped}
    for edge in edges:
        # change edge endpoints if they were matched with a protein in the database
        edge['from'] = nodes_mapped_dict[edge['from']][node_name_attribute] if edge['from'] in nodes_mapped_dict else \
            edge['from']
        edge['to'] = nodes_mapped_dict[edge['to']][node_name_attribute] if edge['to'] in nodes_mapped_dict else edge[
            'to']
    if 'autofill_edges' in parameters['config'] and parameters['config']['autofill_edges']:
        proteins = set(map(lambda n:n[node_name_attribute][1:],filter(lambda n: node_name_attribute in n,parameters['input_network']['nodes'])))
        dataset = 'STRING' if 'interaction_protein_protein' not in parameters['config'] else parameters['config']['interaction_protein_protein']
        dataset_object = models.PPIDataset.objects.filter(name__iexact=dataset).last()
        interaction_objects = models.ProteinProteinInteraction.objects.filter(
            Q(ppi_dataset=dataset_object) & Q(from_protein__in=proteins) & Q(to_protein__in=proteins))
        auto_edges = list(map(lambda n: {"from": f'p{n.from_protein_id}', "to":f'p{n.to_protein_id}'} ,interaction_objects))
        edges.extend(auto_edges)
    result['network']['edges'].extend(edges)

    if not view:
        return Response(result)
    else:
        if view == 'proteins':
            if fmt == 'csv':
                items = []
                for i in proteins:
                    new_i = {
                        'uniprot_ac': i['uniprot_ac'],
                        'gene': i['symbol'],
                        'name': i['protein_name'],
                        'ensg': i['ensg'],
                        'seed': is_seed[i[node_name_attribute]],
                    }
                    if i.get('score'):
                        new_i['score'] = i['score']
                    items.append(new_i)
            else:
                items = proteins
        elif view == 'drugs':
            if fmt == 'csv':
                items = [i for i in drugs]
            else:
                items = drugs
        else:
            return Response({})

        if not fmt or fmt == 'json':
            return Response(items)
        elif fmt == 'csv':
            if len(items) != 0:
                keys = items[0].keys()
            else:
                keys = []
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename="{task.id}_{view}.csv"'
            dict_writer = csv.DictWriter(response, keys)
            dict_writer.writeheader()
            dict_writer.writerows(items)
            return response
        else:
            return Response({})


@api_view(['POST'])
def graph_export(request) -> Response:
    """
    Recieve whole graph data and write it to graphml file. Return the
    file ready to download.
    """
    nodes = request.data.get('nodes', [])
    edges = request.data.get('edges', [])
    fmt = request.data.get('fmt', 'graphml')
    G = nx.Graph()

    for node in nodes:
        # drugstone_id is not interesting outside of drugstone
        # try:
        #     del node['drugstone_id']
        # except KeyError:
        #     pass
        # networkx does not support datatypes such as lists or dicts
        for key in list(node.keys()):
            if isinstance(node[key], list) or isinstance(node[key], dict):
                node[key] = json.dumps(node[key])
            elif node[key] is None:
                # networkx has difficulties with None when writing graphml
                node[key] = ''
        try:
            node_name = node['label']
        except KeyError:
            node_name = node['drugstone_id']
        G.add_node(node_name, **node)

    for e in edges:
        # networkx does not support datatypes such as lists or dicts
        for key in e:
            if isinstance(e[key], list) or isinstance(e[key], dict):
                e[key] = json.dumps(e[key])
            elif e[key] is None:
                e[key] = ''
        u_of_edgece = e.pop('from')
        v_of_edge = e.pop('to')
        G.add_edge(u_of_edgece, v_of_edge, **e)

    if fmt == 'graphml':
        data = nx.generate_graphml(G)
        response = HttpResponse(data, content_type='application/xml')
    elif fmt == 'json':
        data = json.dumps(nx.readwrite.json_graph.node_link_data(G))
        response = HttpResponse(data, content_type='application/json')
    elif fmt == 'csv':
        data = pd.DataFrame(nx.to_numpy_array(G), columns=G.nodes(), index=G.nodes())
        response = HttpResponse(data.to_csv(), content_type='text/csv')

    response['content-disposition'] = f'attachment; filename="{int(time.time())}_network.{fmt}"'
    return response


@api_view(['POST'])
def adjacent_disorders(request) -> Response:
    """Find all adjacent disorders to a list of proteins.

       Args:
           request (django.request): Request object with keys "proteins" and "pdi_dataset"

       Returns:
           Response: With lists "pdis" (protein-drug-intersions) and "disorders"
       """
    data = request.data
    if 'proteins' in data:
        drugstone_ids = data.get('proteins', [])
        pdi_dataset = PDisDataset.objects.filter(name__iexact=data.get('dataset','DisGeNET')).last()
        # find adjacent drugs by looking at drug-protein edges
        pdis_objects = ProteinDisorderAssociation.objects.filter(protein__id__in=drugstone_ids,
                                                                 pdis_dataset=pdi_dataset)
        disorders = {e.disorder for e in pdis_objects}
        # serialize
        edges = ProteinDisorderAssociationSerializer(many=True).to_representation(pdis_objects)
        disorders = DisorderSerializer(many=True).to_representation(disorders)
    elif 'drugs' in data:
        drugstone_ids = data.get('drugs', [])
        drdi_dataset = DrDiDataset.objects.filter(name__iexact=data.get('dataset','DrugBank')).last()
        # find adjacent drugs by looking at drug-protein edges
        drdi_objects = DrugDisorderIndication.objects.filter(drug__id__in=drugstone_ids,
                                                             drdi_dataset=drdi_dataset)
        disorders = {e.disorder for e in drdi_objects}
        # serialize
        edges = DrugDisorderIndicationSerializer(many=True).to_representation(drdi_objects)
        disorders = DisorderSerializer(many=True).to_representation(disorders)
    return Response({
        'edges': edges,
        'disorders': disorders,
    })


@api_view(['POST'])
def adjacent_drugs(request) -> Response:
    """Find all adjacent drugs to a list of proteins.

    Args:
        request (django.request): Request object with keys "proteins" and "pdi_dataset"

    Returns:
        Response: With lists "pdis" (protein-drug-intersions) and "drugs"
    """
    data = request.data
    drugstone_ids = data.get('proteins', [])
    pdi_dataset = PDIDataset.objects.filter(name__iexact=data.get('pdi_dataset','NeDRex')).last()
    # find adjacent drugs by looking at drug-protein edges
    pdi_objects = ProteinDrugInteraction.objects.filter(protein__id__in=drugstone_ids, pdi_dataset=pdi_dataset)
    drugs = {e.drug for e in pdi_objects}
    # serialize
    pdis = ProteinDrugInteractionSerializer(many=True).to_representation(pdi_objects)
    drugs = DrugSerializer(many=True).to_representation(drugs)
    return Response({
        'pdis': pdis,
        'drugs': drugs,
    })


@api_view(['POST'])
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

    return Response({
        'details': details,
        'notFound': not_found,
    })


@api_view(['POST'])
def query_tissue_proteins(request) -> Response:
    threshold = request.data['threshold']
    tissue_id = request.data['tissue_id']
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
        tissue = Tissue.objects.get(id=request.query_params.get('tissue'))

        if request.query_params.get('proteins'):
            ids = json.loads(request.query_params.get('proteins'))
            proteins = list(Protein.objects.filter(id__in=ids).all())
        elif request.query_params.get('token'):
            proteins = []
            task = Task.objects.get(token=request.query_params['token'])
            result = task_result(task)
            network = result['network']
            node_attributes = result.get('node_attributes')
            if not node_attributes:
                node_attributes = {}
            node_types = node_attributes.get('node_types')
            if not node_types:
                node_types = {}
            parameters = json.loads(task.parameters)
            seeds = parameters['seeds']
            nodes = network['nodes']
            for node in nodes + seeds:
                node_type = node_types.get(node)
                details = None
                # if not node_type:
                #     print('we should not see this 3')
                #     node_type, details = infer_node_type_and_details(node)
                if node_type == 'protein':
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
                expression_level = ExpressionLevel.objects.get(protein=protein, tissue=tissue)
                pt_expressions[
                    ProteinSerializer().to_representation(protein)['drugstone_id']] = expression_level.expression_level
            except ExpressionLevel.DoesNotExist:
                pt_expressions[ProteinSerializer().to_representation(protein)['drugstone_id']] = None

        return Response(pt_expressions)
