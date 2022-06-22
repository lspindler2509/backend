# from collections import defaultdict
#
#
# def import_proteins():
#     import python_nedrex as nedrex
#     from python_nedrex.core import get_nodes, get_api_key, get_edges
#     from models import Protein
#
#     def iter_node_collection(coll_name, eval):
#         offset = 0
#         limit = 10000
#         while True:
#             result = get_nodes(coll_name, offset=offset, limit=limit)
#             if not result:
#                 return
#             for node in result:
#                 eval(node)
#             offset += limit
#
#     def iter_edge_collection(coll_name, eval):
#         offset = 0
#         limit = 10000
#         while True:
#             result = get_edges(coll_name, offset=offset, limit=limit)
#             if not result:
#                 return
#             for edge in result:
#                 eval(edge)
#             offset += limit
#
#     def add_protein(node):
#         global proteins
#         id = node['primaryDomainId']
#         proteins[id] = Protein(uniprot_code=id.split('.')[1], gene=node['geneName'])
#
#     def add_edges(edge):
#         global proteins
#         id = edge['sourceDomainId']
#         protein = proteins[id]
#         protein.entrez = edge['targetDomainId'].split('.')[1]
#         global gene_to_prots
#         gene_to_prots[edge['targetDomainId']].add(id)
#
#     def add_genes(node):
#         global proteins
#         global gene_to_prots
#         id = node['primaryDomainId']
#         for prot_id in gene_to_prots[id]:
#             protein = proteins[prot_id]
#             try:
#                 protein.protein_name = node['synonyms'][0]
#             except:
#                 pass
#
#     nedrex.config.set_url_base("http://82.148.225.92:8123/")
#     api_key = get_api_key(accept_eula=True)
#     nedrex.config.set_api_key(api_key)
#
#     proteins = dict()
#     gene_to_prots = defaultdict(lambda: set())
#
#     print('Importing Proteins')
#     iter_node_collection('protein', add_protein)
#     print('Importing Protein-Gene mapping')
#     iter_edge_collection('protein_encoded_by_gene', add_edges)
#     print('Mapping Gene information')
#     iter_node_collection('gene', add_genes)
#     Protein.objects.bulk_create(proteins.values())
