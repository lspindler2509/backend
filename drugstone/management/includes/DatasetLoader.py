from drugstone import models
from python_nedrex.static import get_metadata

def get_ppi_string():
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='STRING',
        link='https://string-db.org/',
        version='11.0'
    )
    return dataset

def get_ppi_apid():
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='APID',
        link='http://cicblade.dep.usal.es:8080/APID/',
        version='January 2019'
    )
    return dataset

def get_ppi_biogrid():
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='BioGRID',
        link='https://thebiogrid.org/',
        version='4.0'
    )
    return dataset

def get_drug_target_nedrex(url):
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='NeDRex',
        link=url,
        version=get_metadata()['version'],
    )
    return dataset

def get_ppi_nedrex(url):
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='NeDRex',
        link=url,
        version=get_metadata()['version'],
    )
    return dataset

def get_protein_disorder_nedrex(url):
    dataset, _ = models.PDisDataset.objects.get_or_create(
        name='NeDRex',
        link=url,
        version=get_metadata()['version'],
    )
    return dataset

def get_drug_disorder_nedrex(url):
    dataset, _ = models.DrDiDataset.objects.get_or_create(
        name='NeDRex',
        link=url,
        version=get_metadata()['version'],
    )
    return dataset

def get_drug_target_chembl():
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='ChEMBL',
        link='https://www.ebi.ac.uk/chembl/',
        version='27',
    )
    return dataset

def get_drug_target_dgidb():
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='DGIdb',
        link='https://www.dgidb.org/',
        version='4.2.0'
    )
    return dataset

def get_drug_target_drugbank():
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='DrugBank',
        link='https://go.drugbank.com/',
        version='5.1.7'
    )
    return dataset

def get_disorder_protein_disgenet():
    dataset, _ = models.PDisDataset.objects.get_or_create(
        name='DisGeNET',
        link='https://www.disgenet.org/home/',
        version='6.0',
    )
    return dataset


def get_drug_disorder_drugbank():
    dataset, _ = models.DrDiDataset.objects.get_or_create(
        name='DrugBank',
        link='https://go.drugbank.com/',
        version='5.1.8',
    )
    return dataset
