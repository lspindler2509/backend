from collections import defaultdict

import python_nedrex as nedrex
from python_nedrex.core import get_nodes, get_edges, get_api_key

from drugstone import models
from drugstone.management.includes.NodeCache import NodeCache


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
            c.append(new_list[id])
        elif new_list[id] != old_list[id]:
            old_list[id].update(new_list[id])
            u.append(old_list[id])
    return u, c


def format_list(l):
    if l is not None and len(l) > 0:
        s = str(l)[1:]
        return s[:len(s) - 1].replace("'", "")
    return ""


def to_id(string):
    idx = string.index('.')
    return string[idx + 1:]


class NedrexImporter:
    cache: NodeCache = None

    def __init__(self, base_url, cache: NodeCache):
        self.cache = cache
        nedrex.config.set_url_base(base_url)
        api_key = get_api_key(accept_eula=True)
        nedrex.config.set_api_key(api_key)

    def import_proteins(self, update: bool):
        proteins = dict()
        gene_to_prots = defaultdict(lambda: set())

        if update:
            self.cache.init_proteins()

        def add_protein(node):
            id = to_id(node['primaryDomainId'])
            name = node['geneName']
            if len(node['synonyms']) > 0:
                name = node['synonyms'][0]
                if '{' in name:
                    idx = name.index('{')
                    if idx > 0:
                        name = name[:idx - 1]
            proteins[id] = models.Protein(uniprot_code=id, protein_name=name, gene=node['geneName'])

        def add_edges(edge):
            id = to_id(edge['sourceDomainId'])
            protein = proteins[id]
            protein.entrez = to_id(edge['targetDomainId'])
            gene_to_prots[protein.entrez].add(id)

        def add_genes(node):
            id = to_id(node['primaryDomainId'])
            for prot_id in gene_to_prots[id]:
                protein = proteins[prot_id]
                try:
                    protein.protein_name = node['synonyms'][0]
                except:
                    pass

        iter_node_collection('protein', add_protein)
        iter_edge_collection('protein_encoded_by_gene', add_edges)

        with_entrez = dict()
        for ids in gene_to_prots.values():
            for id in ids:
                with_entrez[id] = proteins[id]
        proteins = with_entrez

        iter_node_collection('gene', add_genes)

        if update:
            (updates, creates) = identify_updates(proteins, self.cache.proteins)
            for u in updates:
                u.save()
            models.Protein.objects.bulk_create(creates)
            for protein in creates:
                self.cache.proteins[protein.uniprot_code] = protein
                self.cache.protein_updates.add(protein.uniprot_code)
            return len(creates)
        else:
            models.Protein.objects.bulk_create(proteins.values())
            self.cache.proteins = proteins
        return len(self.cache.proteins)

    def import_drugs(self, update):
        drugs = dict()
        if update:
            self.cache.init_drugs()

        def add_drug(node):
            id = to_id(node['primaryDomainId'])
            drugs[id] = models.Drug(drug_id=id, name=node['displayName'], status=format_list(node['drugGroups']))

        iter_node_collection('drug', add_drug)

        if update:
            (updates, creates) = identify_updates(drugs, self.cache.drugs)
            for u in updates:
                u.save()

            models.Drug.objects.bulk_create(creates)
            for drug in creates:
                self.cache.drug_updates.add(drug.drug_id)
                self.cache.drugs[drug.drug_id] = drug
            return len(creates)
        else:
            models.Drug.objects.bulk_create(drugs.values())
            self.cache.drugs = drugs

        return len(self.cache.drugs)

    def import_disorders(self, update):
        disorders = dict()
        if update:
            self.cache.init_disorders()

        def add_disorder(node):
            id = to_id(node['primaryDomainId'])
            disorders[id] = models.Disorder(mondo_id=id, label=node['displayName'], icd10=format_list(node['icd10']))

        iter_node_collection('disorder', add_disorder)

        if update:
            (updates, creates) = identify_updates(disorders, self.cache.disorders)
            for u in updates:
                u.save()
            models.Disorder.objects.bulk_create(creates)
            for disorder in creates:
                self.cache.disorder_updates.add(disorder.mondo_id)
                self.cache.disorders[disorder.mondo_id] = disorder
            return len(creates)
        else:
            models.Disorder.objects.bulk_create(disorders.values())
            self.cache.disorders = disorders

        return len(self.cache.disorders)

    def import_drug_target_interactions(self, dataset, update):
        self.cache.init_drugs()
        self.cache.init_proteins()

        bulk = set()

        def add_dpi(edge):
            try:
                drug = self.cache.get_drug_by_drugbank(to_id(edge['sourceDomainId']))
                protein = self.cache.get_protein_by_uniprot(to_id(edge['targetDomainId']))
                if not update or (self.cache.is_new_drug(drug) or self.cache.is_new_protein(protein)):
                    bulk.add(models.ProteinDrugInteraction(pdi_dataset=dataset, drug=drug, protein=protein))
            except KeyError:
                pass

        iter_edge_collection('drug_has_target', add_dpi)
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def import_protein_protein_interactions(self, dataset, update):
        self.cache.init_proteins()

        bulk = list()

        def iter_ppi(eval):
            from python_nedrex import ppi
            offset = 0
            limit = 10000
            while True:
                result = ppi.ppis({"exp"}, skip=offset, limit=limit)
                if not result:
                    return
                for edge in result:
                    eval(edge)
                offset += limit

        def add_ppi(edge):
            try:
                protein1 = self.cache.get_protein_by_uniprot(to_id(edge['memberOne']))
                protein2 = self.cache.get_protein_by_uniprot(to_id(edge['memberTwo']))
                if not update or (self.cache.is_new_protein(protein1) or self.cache.is_new_protein(protein2)):
                    bulk.append(models.ProteinProteinInteraction(ppi_dataset=dataset, from_protein=protein1,
                                                                 to_protein=protein2))
            except KeyError:
                pass

        iter_ppi(add_ppi)
        models.ProteinProteinInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def import_protein_disorder_associations(self, dataset, update):
        self.cache.init_disorders()
        self.cache.init_proteins()

        bulk = set()

        def add_pdis(edge):
            try:
                disorder = self.cache.get_disorder_by_mondo(to_id(edge['targetDomainId']))
                for protein in self.cache.get_proteins_by_entrez(to_id(edge['sourceDomainId'])):
                    if not update or (self.cache.is_new_disease(disorder) or self.cache.is_new_protein(protein)):
                        bulk.add(models.ProteinDisorderAssociation(pdis_dataset=dataset, protein=protein,
                                                                   disorder=disorder, score=edge['score']))
            except KeyError:
                pass

        iter_edge_collection('gene_associated_with_disorder', add_pdis)
        models.ProteinDisorderAssociation.objects.bulk_create(bulk)
        return len(bulk)

    def import_drug_disorder_indications(self, dataset, update):
        self.cache.init_disorders()
        self.cache.init_drugs()

        bulk = set()

        def add_drdis(edge):
            try:
                drug = self.cache.get_drug_by_drugbank(to_id(edge['sourceDomainId']))
                disorder = self.cache.get_disorder_by_mondo(to_id(edge['targetDomainId']))
                if not update or (self.cache.is_new_drug(drug) or self.cache.is_new_disease(disorder)):
                    bulk.add(models.DrugDisorderIndication(drdi_dataset=dataset, drug=drug, disorder=disorder))
            except KeyError:
                pass

        iter_edge_collection('drug_has_indication', add_drdis)
        models.DrugDisorderIndication.objects.bulk_create(bulk)
        return len(bulk)
