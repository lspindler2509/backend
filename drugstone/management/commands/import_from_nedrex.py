from collections import defaultdict

import python_nedrex as nedrex
from python_nedrex.core import get_nodes, get_edges, get_api_key

from drugstone import models


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


def identify_updates(new_list, old_list):
    u = list()
    c = list()
    for id in new_list:
        if id not in old_list:
            c.append(id)
        elif new_list[id] != old_list[id]:
            old_list[id].update(new_list[id])
            u.append(old_list[id])
    return u, c


def format_list(l):
    if l is not None and len(l) > 0:
        s = str(l)[1:]
        return s[:len(s) - 1]
    return ""


class nedrex_importer:
    proteins = dict()
    entrez_to_uniprot = dict()
    gene_name_to_uniprot = defaultdict(lambda: set())
    disorders = dict()
    drugs = dict()

    def __init__(self, base_url):
        nedrex.config.set_url_base(base_url)
        api_key = get_api_key(accept_eula=True)
        nedrex.config.set_api_key(api_key)

    def init_proteins(self):
        if len(self.proteins) == 0:
            print("Generating protein maps...")
            for protein in models.Protein.objects.all():
                self.proteins[protein.entrez] = protein
                self.entrez_to_uniprot[protein.entrez] = protein.uniprot_code
                self.gene_name_to_uniprot[protein.gene].add(protein.uniprot_code)

    def init_drugs(self):
        if len(self.drugs) == 0:
            print("Generating drug map...")
            for drug in models.Drug.objects.all():
                self.drugs[drug.drug_id] = drug

    def init_disorders(self):
        if len(self.disorders) == 0:
            print("Generating disorder map...")
            for disorder in models.Disorder.objects.all():
                self.disorders[disorder.mondo_id] = disorder

    def import_proteins(self, update: bool):
        proteins = dict()
        gene_to_prots = defaultdict(lambda: set())

        if update:
            self.init_proteins()

        def add_protein(node):
            print(node)
            id = node['primaryDomainId'].split('.')[1]
            name = node['geneName']
            if len(node['synonyms']) > 0:
                name = node['synonyms'][0]
                idx = name.index('{')
                if idx > 0:
                    name = name[idx - 1:]
            proteins[id] = models.Protein(uniprot_code=id, name=name, gene=node['geneName'])

        def add_edges(edge):
            id = edge['sourceDomainId'].split('.')[1]
            protein = proteins[id]
            protein.entrez = edge['targetDomainId'].split('.')[1]
            gene_to_prots[edge['targetDomainId']].add(id)

        def add_genes(node):
            id = node['primaryDomainId'].split('.')[1]
            for prot_id in gene_to_prots[id]:
                protein = proteins[prot_id]
                try:
                    protein.protein_name = node['synonyms'][0]
                except:
                    pass

        iter_node_collection('protein', add_protein)
        iter_edge_collection('protein_encoded_by_gene', add_edges)
        iter_node_collection('gene', add_genes)
        # TODO test updating ideas
        if update:
            (updates, creates) = identify_updates(proteins, self.proteins)
            models.Protein.objects.bulk_update(updates)
            models.Protein.objects.bulk_create(creates)
            for protein in creates:
                self.proteins[protein.uniprot_code] = protein
        else:
            models.Protein.objects.bulk_create(self.proteins.values())
            self.proteins = proteins
        return len(self.proteins)

    def import_drugs(self, update):
        drugs = dict()
        if update:
            self.init_drugs()

        def add_drug(node):
            id = node['primaryDomainId'].split('.')[1]
            drugs[id] = models.Drug(drug_id=id, name=node['displayName'], status=format_list(node['drugGroups']))

        iter_node_collection('drug', add_drug)

        # TODO test updating ideas
        if update:
            (updates, creates) = identify_updates(drugs, self.drugs)
            models.Drug.objects.bulk_update(updates)
            models.Drug.objects.bulk_create(creates)
            for drug in creates:
                self.drugs[drug.drug_id] = drug
        else:
            models.Drug.objects.bulk_create(self.drugs.values())
            self.drugs = drugs

        self.drugs = drugs
        return len(self.drugs)

    def import_disorders(self, update):
        disorders = dict()
        if update:
            self.init_disorders()

        def add_disorder(node):
            id = node['primaryDomainId'].split('.')[1]
            self.disorders[id] = models.Disorder(mondo_id=id, label=node['displayName'], icd10=format_list(node['icd10']))

        iter_node_collection('disorder', add_disorder)

        # TODO test updating ideas
        if update:
            (updates, creates) = identify_updates(disorders, self.disorders)
            models.Disorder.objects.bulk_update(updates)
            models.Disorder.objects.bulk_create(creates)
            for disorder in creates:
                self.disorders[disorder.uniprot_code] = disorder
        else:
            models.Disorder.objects.bulk_create(self.disorders.values())
            self.disorders = disorders

        self.disorders = disorders
        return len(self.disorders)
