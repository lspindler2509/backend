import copy
from collections import defaultdict
import json
from typing import List, Tuple, Set, OrderedDict
from functools import reduce
from django.db.models import Q
from drugstone.models import Protein, EnsemblGene, Task
from drugstone.serializers import ProteinSerializer


MAP_ID_SPACE_COMPACT_TO_DRUGSTONE = {
    'symbol:': 'symbol',
    'uniprot:': 'uniprot',
    'ensg:': 'ensg',
    'ncbigene:': 'entrez',
    'ensembl:': 'ensg',
    'entrez:': 'entrez'
}


def query_proteins_by_identifier(node_ids: Set[str], identifier: str) -> Tuple[List[dict], str]:
    """Queries the django database Protein table given a list of identifiers (node_ids) and a identifier name
    (identifier).
    The identifier name represents any protein attribute, e.g. uniprot or symbol.
    The identifier names vary from the Protein table names since they are the strings which are set by the user
    in the frontend, for readability they were changes from the original backend attributes.

    Args:
        node_ids (list): List of protein or gene identifiers. Note: Do not mix identifiers.
        identifier (str): Can be one of "symbol", "ensg", "uniprot"

    Returns:
        Tuple[List[dict], str]:
            Returns list of serialized protein entries for all matched IDs
            Returns name of backend attribute of Protein table
    """
    # query protein table
    if (len(node_ids) == 0):
        return list(), identifier
    if identifier == 'symbol':
        protein_attribute = 'symbol'
        q_list = map(lambda n: Q(gene__iexact=n), node_ids)
    elif identifier == 'uniprot':
        protein_attribute = 'uniprot'
        q_list = map(lambda n: Q(uniprot_code__iexact=n), node_ids)
    elif identifier == 'ensg' or identifier == 'ensembl':
        protein_attribute = 'ensg'
        dr_ids = map(lambda n: n.protein_id, EnsemblGene.objects.filter(
            reduce(lambda a, b: a | b, map(lambda n: Q(name__iexact=n), list(node_ids)))))
        q_list = map(lambda n: Q(id=n), dr_ids)
    elif identifier == 'entrez' or identifier == 'ncbigene':
        protein_attribute = 'entrez'
        q_list = map(lambda n: Q(entrez=n), node_ids)
    if not node_ids:
        # node_ids is an empty list
        return [], protein_attribute
    q_list = reduce(lambda a, b: a | b, q_list)
    node_objects = Protein.objects.filter(q_list)
    
    cp_to_node_ids = {}
    for node in node_objects:
        cellular_components = node.cellular_components.all()
        components = []
        for cp in cellular_components:
            cp_string = cp.go_code + ":" + cp.display_name
            components.append(cp_string)
        
        node_id = ''
        if protein_attribute == 'symbol':
            node_id = node.gene
        elif protein_attribute == 'uniprot':
            node_id = node.uniprot_code
        elif protein_attribute == 'entrez':
            node_id = node.entrez
        cp_to_node_ids[node_id] = components
        
    nodes = list()
    node_map = defaultdict(list)
    if protein_attribute == 'ensg':
        for node in ProteinSerializer(many=True).to_representation(node_objects):
            for ensembl_id in node.get(protein_attribute):
                if ensembl_id.upper() in node_ids:
                    node = copy.copy(node)
                    node[identifier] = ensembl_id
                    node_map[ensembl_id].append(node)
    else:
        for node in ProteinSerializer(many=True).to_representation(node_objects):
            id_node = node.get(protein_attribute)
            if id_node in cp_to_node_ids:
                node["cellular_component"] = cp_to_node_ids[id_node]
            else:
                node["cellular_component"] = []
            node_map[node.get(protein_attribute)].append(node)
    for node_id, entries in node_map.items():
        nodes.append(aggregate_nodes(entries))
    return nodes, protein_attribute


def get_protein_ids(id_space, proteins):
    if (id_space == 'uniprot'):
        return {p['uniprot'] for p in proteins}
    if (id_space == 'ensg' or id_space == 'ensembl'):
        return {p['ensg'] for p in proteins}
    if (id_space == 'symbol'):
        return {p['symbol'] for p in proteins}
    if (id_space == 'entrez' or id_space == 'ncbigene'):
        return {p['entrez'] for p in proteins}
    return set()


def clean_proteins_from_compact_notation(node_ids: Set[str], identifier: str) -> List[str]:
    """Queries the django database Protein table given a list of identifiers (node_ids) and a identifier name
    (identifier).
    The identifier name represents any protein attribute, e.g. uniprot or symbol.
    The identifier names vary from the Protein table names since they are the strings which are set by the user
    in the frontend, for readability they were changes from the original backend attributes.

    Args:
        node_ids (list): List of protein or gene identifiers. Note: Do not mix identifiers.
        identifier (str): Can be one of "symbol", "ensg", "uniprot"

    Returns:
        Tuple[List[dict], str]:
            Returns list of serialized protein entries for all matched IDs
            Returns name of backend attribute of Protein table
    """
    # query protein table
    if len(node_ids) == 0:
        return list()

    symbol_set, ensg_set, uniprot_set, entrez_set = set(), set(), set(), set()

    id_map = {
        'symbol:': symbol_set,
        'uniprot:': uniprot_set,
        'ensg:': ensg_set,
        'ncbigene:': entrez_set,
        'ensembl:': ensg_set,
        'entrez:': entrez_set
    }
    clean_ids = set()
    for node_id in node_ids:
        added = False
        for id_space in id_map.keys():
            if node_id.startswith(id_space):
                id_map[id_space].add(node_id[len(id_space):].upper())
                added = True
                break
        if not added:
            clean_ids.add(node_id)

    for id_space, ids in id_map.items():
        if len(ids) == 0:
            continue
        if id_space == 'symbol:':
            q_list = map(lambda n: Q(gene__iexact=n), ids)
        elif id_space == 'uniprot:':
            q_list = map(lambda n: Q(uniprot_code__iexact=n), ids)
        elif id_space == 'ensg:':
            ensembls = EnsemblGene.objects.filter(reduce(lambda a, b: a | b, map(lambda n: Q(name__iexact=n), ids)))
            if len(ensembls) == 0:
                continue
            dr_ids = map(lambda n: n.protein_id, ensembls)
            q_list = map(lambda n: Q(id=n), dr_ids)
        elif id_space == 'entrez:':
            q_list = map(lambda n: Q(entrez=n), ids)
        else:
            continue
        q_list = reduce(lambda a, b: a | b, q_list)
        proteins = ProteinSerializer(many=True).to_representation(Protein.objects.filter(q_list))
        # if protein could not be mapped
        clean_ids_temp = get_protein_ids(identifier, proteins)
        if '' in clean_ids_temp:
            clean_ids_temp.remove('')
            # at least one protein could not be found in id space, use original id as placeholder
            ids_placeholder = {p[MAP_ID_SPACE_COMPACT_TO_DRUGSTONE[id_space]] for p in proteins if p[identifier] == ''}
            clean_ids_temp |= ids_placeholder
        clean_ids |= clean_ids_temp

    return list(clean_ids)


def aggregate_nodes(nodes: List[OrderedDict]):
    node = defaultdict(set)
    for n in nodes:
        for key, value in n.items():
            if isinstance(value, list):
                for e in value:
                    if e is not None and len(e) > 0:
                        node[key].add(e)
            elif value is not None and len(value) > 0:
                node[key].add(value)
    return {k: list(v) for k, v in node.items()}

def update_result(result, token: str):
    task = Task.objects.get(token=token)
    task.result = json.dumps(result)
    task.save()
    
