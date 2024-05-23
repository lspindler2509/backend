import json
from collections import defaultdict

import nedrex
from nedrex.core import get_nodes, get_edges, get_api_key

from drugstone import models
from drugstone.management.includes.NodeCache import NodeCache
from drugstone.management.includes import DatasetLoader
from drugstone.models import PPIDataset


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
    url: str = ''
    licenced_url: str = ''
    unlicenced_url: str = ''
    licenced_on: bool = True
    api_key: str = None

    def __init__(self, base_url_licenced, base_url_unlicenced, cache: NodeCache):
        self.cache = cache
        self.licenced_url = base_url_licenced
        self.unlicenced_url = base_url_unlicenced
        self.set_licenced(False)

    def get_api_key(self):
        if self.api_key is None:
            self.api_key = get_api_key(accept_eula=True)
        return self.api_key

    def set_licenced(self, on):
        if on == self.licenced_on:
            return

        self.url = self.licenced_url if on else self.unlicenced_url
        nedrex.config.set_url_base(self.url)

        if on:
            nedrex.config.set_api_key(self.get_api_key())

        self.licenced_on = on

    def import_proteins(self, update: bool):
        self.set_licenced(False)
        proteins = dict()
        gene_to_prots = defaultdict(lambda: set())

        if update:
            self.cache.init_proteins()

        def format_prot_name(name):
            if '{' in name:
                idx1 = name.index('{')
                adjusted_name = name[:idx1 - 1].strip() if idx1 > 0 else ''
                if '=' in adjusted_name:
                    idx2 = adjusted_name.index('=')
                    return adjusted_name[idx2 + 1:].strip()
                return adjusted_name
            return name

        def add_protein(node):
            id = to_id(node['primaryDomainId'])
            name = format_prot_name(node['geneName'])
            gene = name

            if len(node['synonyms']) > 0:
                name = format_prot_name(node['synonyms'][0])
            proteins[id] = models.Protein(uniprot_code=id, protein_name=name, gene=gene)

        def add_edges(edge):
            try:
                id = to_id(edge['sourceDomainId'])
                protein = proteins[id]
                protein.entrez = to_id(edge['targetDomainId'])
                gene_to_prots[protein.entrez].add(id)
            except:
                print(f'Edge could not be saved: {edge["sourceDomainId"]} - {edge["targetDomainId"]}')

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
                self.cache.proteins[u.uniprot_code] = u
            models.Protein.objects.bulk_create(creates)
            for protein in creates:
                self.cache.proteins[protein.uniprot_code] = protein
                self.cache.protein_updates.add(protein.uniprot_code)
            self.cache.init_protein_maps()
            return len(creates)
        else:
            models.Protein.objects.bulk_create(proteins.values())
            self.cache.proteins = proteins
        return len(self.cache.proteins)

    def import_drugs(self, update):
        self.set_licenced(False)

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
        licenced = dataset.licenced
        self.set_licenced(licenced)

        self.cache.init_drugs()
        self.cache.init_proteins()

        bulk = set()
        delete = set()
        existing = dict()
        if update:
            for edge in models.ProteinDrugInteraction.objects.filter(pdi_dataset=dataset):
                existing[edge.__hash__()] = edge

        source_datasets = dict()
        source_is_licenced = dict()

        def get_dataset(source):
            if source not in source_datasets:
                source_datasets[source] = DatasetLoader.get_pdi_nedrex_dataset(self.url, licenced, source)
            return source_datasets[source]

        def is_licenced(source):
            if source not in source_is_licenced:
                source_is_licenced[source] = DatasetLoader.is_licenced_pdi_source(source)
            return source_is_licenced[source]

        def add_dpi(edge):
            try:
                drug = self.cache.get_drug_by_drugbank(to_id(edge['sourceDomainId']))
                protein = self.cache.get_protein_by_uniprot(to_id(edge['targetDomainId']))
                actions = json.dumps(edge['actions'])
                e = models.ProteinDrugInteraction(pdi_dataset=dataset, drug=drug, protein=protein, actions=actions)
                if update and e.__hash__() in existing:
                    if existing[e.__hash__()] != e:
                        delete.add(existing[e.__hash__()])
                        del existing[e.__hash__()]
                    else:
                        return
                if not update or e.__hash__() not in existing:
                    bulk.add(e)
                    for source in edge['dataSources']:
                        if licenced:
                            if not is_licenced(source):
                                continue
                        bulk.add(models.ProteinDrugInteraction(pdi_dataset=get_dataset(source), drug=drug,
                                                               protein=protein))
            except KeyError:
                pass

        iter_edge_collection('drug_has_target', add_dpi)
        for d in delete:
            d.delete()
        
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def import_protein_protein_interactions(self, dataset: PPIDataset, update):
        licenced = dataset.licenced
        self.set_licenced(licenced)

        self.cache.init_proteins()

        bulk = list()
        existing = set()
        if update:
            for edge in models.ProteinProteinInteraction.objects.filter(ppi_dataset=dataset):
                existing.add(edge.__hash__())

        source_datasets = dict()
        source_is_licenced = dict()

        def get_dataset(source):
            if source not in source_datasets:
                source_datasets[source] = DatasetLoader.get_ppi_nedrex_dataset(self.url, licenced, source)
            return source_datasets[source]

        def is_licenced(source):
            if source not in source_is_licenced:
                source_is_licenced[source] = DatasetLoader.is_licenced_ppi_source(source)
            return source_is_licenced[source]

        def iter_ppi(eval):
            from nedrex import ppi
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
                e = models.ProteinProteinInteraction(ppi_dataset=dataset, from_protein=protein1, to_protein=protein2)
                if not update or e.__hash__() not in existing:
                    bulk.append(e)
                    for source in edge['dataSources']:
                        if licenced:
                            if not is_licenced(source):
                                continue
                        bulk.append(
                            models.ProteinProteinInteraction(ppi_dataset=get_dataset(source), from_protein=protein1,
                                                             to_protein=protein2))
            except KeyError:
                pass

        iter_ppi(add_ppi)
        models.ProteinProteinInteraction.objects.bulk_create(bulk)
        # new_datasets = [dataset, source_datasets.values()]
        # DatasetLoader.remove_old_ppi_data(new_datasets, licenced)
        return len(bulk)

    def import_protein_disorder_associations(self, dataset, update):
        licenced = dataset.licenced
        self.set_licenced(licenced)

        self.cache.init_disorders()
        self.cache.init_proteins()

        bulk = set()
        existing = set()
        if update:
            for edge in models.ProteinDisorderAssociation.objects.filter(pdis_dataset=dataset):
                existing.add(edge.__hash__())

        source_datasets = dict()
        source_is_licenced = dict()

        def get_dataset(source):
            if source not in source_datasets:
                source_datasets[source] = DatasetLoader.get_pdis_nedrex_dataset(self.url, licenced, source)
            return source_datasets[source]

        def is_licenced(source):
            if source not in source_is_licenced:
                source_is_licenced[source] = DatasetLoader.is_licenced_pdis_source(source)
            return source_is_licenced[source]

        def add_pdis(edge):
            try:
                disorder = self.cache.get_disorder_by_mondo(to_id(edge['targetDomainId']))
                for protein in self.cache.get_proteins_by_entrez(to_id(edge['sourceDomainId'])):
                    e = models.ProteinDisorderAssociation(pdis_dataset=dataset, protein=protein, disorder=disorder,
                                                          score=edge['score'])
                    if not update or e.__hash__() not in existing:
                        bulk.add(e)
                        for source in edge['dataSources']:
                            if licenced:
                                if not is_licenced(source):
                                    continue
                            bulk.add(
                                models.ProteinDisorderAssociation(pdis_dataset=get_dataset(source), protein=protein,
                                                                  disorder=disorder,
                                                                  score=edge['score']))
            except KeyError:
                pass

        iter_edge_collection('gene_associated_with_disorder', add_pdis)
        models.ProteinDisorderAssociation.objects.bulk_create(bulk)
        # new_datasets = [dataset, source_datasets.values()]
        # DatasetLoader.remove_old_pdis_data(new_datasets, licenced)
        return len(bulk)

    def import_drug_disorder_indications(self, dataset, update):
        licenced = dataset.licenced
        self.set_licenced(licenced)

        self.cache.init_disorders()
        self.cache.init_drugs()

        bulk = set()
        existing = set()
        if update:
            for edge in models.DrugDisorderIndication.objects.filter(drdi_dataset=dataset):
                existing.add(edge.__hash__())

        source_datasets = dict()
        source_is_licenced = dict()

        def get_dataset(source):
            if source not in source_datasets:
                source_datasets[source] = DatasetLoader.get_drdi_nedrex_dataset(self.url, licenced, source)
            return source_datasets[source]

        def is_licenced(source):
            if source not in source_is_licenced:
                source_is_licenced[source] = DatasetLoader.is_licenced_drdi_source(source)
            return source_is_licenced[source]

        def add_drdis(edge):
            try:
                drug = self.cache.get_drug_by_drugbank(to_id(edge['sourceDomainId']))
                disorder = self.cache.get_disorder_by_mondo(to_id(edge['targetDomainId']))
                e = models.DrugDisorderIndication(drdi_dataset=dataset, drug=drug, disorder=disorder)
                if not update or e.__hash__() not in existing:
                    bulk.add(e)
                    for source in edge['dataSources']:
                        if licenced:
                            if not is_licenced(source):
                                continue
                        bulk.add(
                            models.DrugDisorderIndication(drdi_dataset=get_dataset(source), drug=drug,
                                                          disorder=disorder))
            except KeyError:
                return

        iter_edge_collection('drug_has_indication', add_drdis)
        models.DrugDisorderIndication.objects.bulk_create(bulk)
        # new_datasets = [dataset, source_datasets.values()]
        # DatasetLoader.remove_old_drdi_data(new_datasets, licenced)
        return len(bulk)
