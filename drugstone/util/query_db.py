from typing import List, Tuple
from functools import reduce
from django.db.models import Q
from drugstone.models import Protein
from drugstone.serializers import ProteinSerializer


def query_proteins_by_identifier(node_ids: List[str], identifier: str) -> Tuple[List[dict], str]:
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
        q_list = map(lambda n: Q(ensg__name__iexact=n), node_ids)

    if not node_ids:
        # node_ids is an empty list
        return [], protein_attribute

    q_list = reduce(lambda a, b: a | b, q_list)
    node_objects = Protein.objects.filter(q_list)
    # serialize
    nodes = ProteinSerializer(many=True).to_representation(node_objects)

    return nodes, protein_attribute
