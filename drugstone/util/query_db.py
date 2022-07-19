from collections import defaultdict
from typing import List, Tuple, Set, OrderedDict
from functools import reduce
from django.db.models import Q
from drugstone.models import Protein, EnsemblGene
from drugstone.serializers import ProteinSerializer


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
    if identifier == 'symbol':
        protein_attribute = 'symbol'
        q_list = map(lambda n: Q(gene__iexact=n), node_ids)
    elif identifier == 'uniprot':
        protein_attribute = 'uniprot_ac'
        q_list = map(lambda n: Q(uniprot_code__iexact=n), node_ids)
    elif identifier == 'ensg':
        protein_attribute = 'ensg'
        node_ids = map(lambda n: n.protein_id, EnsemblGene.objects.filter(
            reduce(lambda a, b: a | b, map(lambda n: Q(name__iexact=n), list(node_ids)))))
        q_list = map(lambda n: Q(id=n), node_ids)
    elif identifier == 'entrez':
        protein_attribute = 'entrez'
        q_list = map(lambda n: Q(entrez=n), node_ids)
    if not node_ids:
        # node_ids is an empty list
        return [], protein_attribute
    q_list = reduce(lambda a, b: a | b, q_list)
    node_objects = Protein.objects.filter(q_list)

    nodes = list()

    node_map = defaultdict(list)

    for node in ProteinSerializer(many=True).to_representation(node_objects):
        node_map[node.get(protein_attribute)].append(node)
    for node_id, entries in node_map.items():
        nodes.append(aggregate_nodes(entries))

    return nodes, protein_attribute


def aggregate_nodes(nodes: List[OrderedDict]):
    node = defaultdict(set)
    for n in nodes:
        for key, value in n.items():
            if isinstance(value,list):
                for e in value:
                    node[key].add(e)
            else:
                node[key].add(value)
    return {k: list(v) for k, v in node.items()}
