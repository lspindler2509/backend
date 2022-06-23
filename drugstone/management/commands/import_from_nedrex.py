from collections import defaultdict


class nedrex_importer:

    proteins = dict()

    def import_proteins(self):
        import python_nedrex as nedrex
        from python_nedrex.core import get_nodes, get_api_key, get_edges
        from drugstone.models import Protein

        gene_to_prots = defaultdict(lambda: set())

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

        def add_protein(node):
            id = node['primaryDomainId']
            self.proteins[id] = Protein(uniprot_code=id.split('.')[1], gene=node['geneName'])

        def add_edges(edge):
            id = edge['sourceDomainId']
            protein = self.proteins[id]
            protein.entrez = edge['targetDomainId'].split('.')[1]
            gene_to_prots[edge['targetDomainId']].add(id)

        def add_genes(node):
            id = node['primaryDomainId']
            for prot_id in gene_to_prots[id]:
                protein = self.proteins[prot_id]
                try:
                    protein.protein_name = node['synonyms'][0]
                except:
                    pass

        nedrex.config.set_url_base("http://82.148.225.92:8123/")
        api_key = get_api_key(accept_eula=True)
        nedrex.config.set_api_key(api_key)



        print('Importing Proteins')
        iter_node_collection('protein', add_protein)
        print('Importing Protein-Gene mapping')
        iter_edge_collection('protein_encoded_by_gene', add_edges)
        print('Mapping Gene information')
        iter_node_collection('gene', add_genes)
        Protein.objects.bulk_create(self.proteins.values())
        return len(self.proteins)
