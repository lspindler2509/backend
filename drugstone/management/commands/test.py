import python_nedrex as nedrex
from python_nedrex.core import get_nodes, get_edges, get_api_key
from python_nedrex.static import get_metadata

def iter_node_collection(coll_name, eval):
    offset = 0
    limit = 10000
    while True:
        result = get_nodes(coll_name, offset=offset, limit=limit)
        if not result:
            return
        for node in result:
            eval(node)
        offset += limit


def iter_edge_collection(coll_name, eval):
    offset = 0
    limit = 10000
    while True:
        result = get_edges(coll_name, offset=offset, limit=limit)
        if not result:
            return
        for edge in result:
            eval(edge)
        offset += limit


def iter_ppi(eval):
    from python_nedrex import ppi
    offset = 0
    limit = 1000
    while True:
        result = ppi.ppis({"exp"},skip = offset, limit=limit)
        if not result:
            return
        for edge in result:
            eval(edge)
        offset += limit

base_url = "http://82.148.225.92:8123/"
nedrex.config.set_url_base(base_url)
api_key = get_api_key(accept_eula=True)
nedrex.config.set_api_key(api_key)
print(f'Nodes: {nedrex.core.get_node_types()}')
print(f'Edges: {nedrex.core.get_edge_types()}')
print(f'{get_metadata()}')


iter_ppi(lambda node: print(node))
# iter_edge_collection("gene_expressed_in_tissue", lambda node: {print(node)})