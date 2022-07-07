from drugstone import models
from python_nedrex.static import get_metadata

ppi_nedrex_datasets = dict()


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


def get_ppi_nedrex_biogrid(url):
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='BioGRID',
        link=url,
        version=get_metadata()['source_databases']['biogrid']['date']
    )
    return dataset


def get_ppi_nedrex_iid(url):
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='IID',
        link=url,
        version=get_metadata()['source_databases']['iid']['date']
    )
    return dataset


def get_ppi_nedrex_intact(url):
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='IntAct',
        link=url,
        version=get_metadata()['source_databases']['intact']['date']
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


def get_dis_prot_nedrex_disgenet(url):
    dataset, _ = models.PDisDataset.objects.get_or_create(
        name='DisGeNET',
        link=url,
        version=get_metadata()['source_databases']['disgenet']['date']
    )
    return dataset


def get_dis_prot_nedrex_omim(url):
    dataset, _ = models.PDisDataset.objects.get_or_create(
        name='OMIM',
        link=url,
        version=get_metadata()['source_databases']['omim']['date']
    )
    return dataset


def get_drdis_nedrex_drugcentral(url):
    dataset, _ = models.DrDiDataset.objects.get_or_create(
        name='Drug Central',
        link=url,
        version=get_metadata()['source_databases']['drug_central']['date']
    )
    return dataset

def get_drdis_nedrex_ctd(url):
    dataset, _ = models.DrDiDataset.objects.get_or_create(
        name='CTD',
        link=url,
        version=get_metadata()['source_databases']['ctd']['date']
    )
    return dataset

def get_pdr_nedrex_drugcentral(url):
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='Drug Central',
        link=url,
        version=get_metadata()['source_databases']['drug_central']['date']
    )
    return dataset

def get_pdr_nedrex_drugbank(url):
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='DrugBank',
        link=url,
        version=get_metadata()['source_databases']['drugbank']['date']
    )
    return dataset

def get_pdr_nedrex_datasets(url):
    return {'DrugBank': get_pdr_nedrex_drugbank(url), 'DrugCentral': get_pdr_nedrex_drugcentral(url)}

def get_drdis_nedrex_datasets(url):
    return {'ctd':get_drdis_nedrex_ctd(url), 'Drug Central':get_drdis_nedrex_drugcentral(url)}

def get_ppi_nedrex_datasets(url):
    return {'biogrid':get_ppi_nedrex_biogrid(url), 'iid':get_ppi_nedrex_iid(url), 'intact':get_ppi_nedrex_intact(url)}

def get_dis_prot_nedrex_datasets(url):
    return {'disgenet': get_dis_prot_nedrex_disgenet(url), 'omim': get_dis_prot_nedrex_omim(url)}