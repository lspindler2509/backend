from requests.exceptions import RetryError

from drugstone import models
from nedrex.static import get_metadata, get_license

LICENSE_FILE = "./data/license.txt"


def get_ppi_string():
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='STRING',
        link='https://string-db.org/',
        version='11.0',
        licenced=False
    )
    return dataset


def get_ppi_apid():
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='APID',
        link='http://cicblade.dep.usal.es:8080/APID/',
        version='January 2019',
        licenced=False
    )
    return dataset


def get_ppi_biogrid():
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='BioGRID',
        link='https://thebiogrid.org/',
        version='4.0',
        licenced=False
    )
    return dataset


def get_nedrex_version():
    version = get_today_version()
    try:
        version = get_metadata()['version']
    except RetryError:
        pass
    return version


def get_nedrex_source_version(source):
    metadata = get_metadata()['source_databases']
    # TODO remove once fixed in nedrex db
    if 'drug_central' in metadata:
        metadata['drugcentral'] = metadata['drug_central']
        
    return metadata[source]['date']


def get_drug_target_nedrex(url, licenced):
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='NeDRex',
        link=url,
        version=get_nedrex_version(),
        licenced=licenced
    )
    return dataset


def get_ppi_nedrex(url, licenced):
    dataset, _ = models.PPIDataset.objects.get_or_create(
        name='NeDRex',
        link=url,
        version=get_nedrex_version(),
        licenced=licenced
    )
    return dataset


def get_protein_disorder_nedrex(url, licenced):
    dataset, _ = models.PDisDataset.objects.get_or_create(
        name='NeDRex',
        link=url,
        version=get_nedrex_version(),
        licenced=licenced
    )
    return dataset


def get_drug_disorder_nedrex(url, licenced):
    dataset, _ = models.DrDiDataset.objects.get_or_create(
        name='NeDRex',
        link=url,
        version=get_nedrex_version(),
        licenced=licenced
    )
    return dataset


def write_license(text):
    with open(LICENSE_FILE, 'w') as fh:
        fh.write(text)


def update_license():
    try:
        license = get_license()
        write_license(license)
        return license
    except RetryError:
        print(f'License could not be retreived.')
        return ""


def import_license():
    try:
        license = ""
        with open(LICENSE_FILE, 'r') as fh:
            for line in fh:
                license += line
        return license
    except FileNotFoundError:
        print(f'No license doc there yet! Make sure to run an update first!')
    return ""


def get_drug_target_chembl():
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='ChEMBL',
        link='https://www.ebi.ac.uk/chembl/',
        version='27',
        licenced=False
    )
    return dataset


def get_drug_target_dgidb():
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='DGIdb',
        link='https://www.dgidb.org/',
        version='4.2.0',
        licenced=False
    )
    return dataset


def get_drug_target_drugbank():
    dataset, _ = models.PDIDataset.objects.get_or_create(
        name='DrugBank',
        link='https://go.drugbank.com/',
        version='5.1.7',
        licenced=True
    )
    return dataset


def get_disorder_protein_disgenet():
    dataset, _ = models.PDisDataset.objects.get_or_create(
        name='DisGeNET',
        link='https://www.disgenet.org/home/',
        version='6.0',
        licenced=False
    )
    return dataset


def get_drug_disorder_drugbank():
    dataset, _ = models.DrDiDataset.objects.get_or_create(
        name='DrugBank',
        link='https://go.drugbank.com/',
        version='5.1.8',
        licenced=False
    )
    return dataset


def get_today_version():
    import datetime
    now = datetime.date.today()
    version = f'{now.year}-{now.month}-{now.day}_temp'
    return version


def get_ppi_nedrex_dataset(url, licenced, source):
    version = get_today_version()
    try:
        version = get_nedrex_source_version(source)
    except RetryError:
        pass

    dataset, _ = models.PPIDataset.objects.get_or_create(
        name=source,
        link=url,
        version=version,
        licenced=licenced
    )
    return dataset


def get_pdi_nedrex_dataset(url, licenced, source):
    version = get_today_version()
    try:
        version = get_nedrex_source_version(source)
    except RetryError:
        pass

    dataset, _ = models.PDIDataset.objects.get_or_create(
        name=source,
        link=url,
        version=version,
        licenced=licenced
    )
    return dataset


def get_pdis_nedrex_dataset(url, licenced, source):
    version = get_today_version()
    try:
        version = get_nedrex_source_version(source)
    except RetryError:
        pass

    dataset, _ = models.PDisDataset.objects.get_or_create(
        name=source,
        link=url,
        version=version,
        licenced=licenced
    )
    return dataset


def get_drdi_nedrex_dataset(url, licenced, source):
    version = get_today_version()
    try:
        version = get_nedrex_source_version(source)
    except RetryError:
        pass

    dataset, _ = models.DrDiDataset.objects.get_or_create(
        name=source,
        link=url,
        version=version,
        licenced=licenced
    )
    return dataset


def is_licenced_ppi_source(source):
    version = get_today_version()
    try:
        version = get_nedrex_source_version(source)
    except RetryError:
        pass

    try:
        models.PPIDataset.objects.get(name=source, version=version, licenced=False).link
    except:
        return True
    return False


def is_licenced_pdi_source(source):
    version = get_today_version()
    try:
        version = get_nedrex_source_version(source)
    except RetryError:
        pass

    try:
        models.PDIDataset.objects.get(name=source, version=version, licenced=False).link
    except:
        return True
    return False


def is_licenced_pdis_source(source):
    version = get_today_version()
    try:
        version = get_nedrex_source_version(source)
    except RetryError:
        pass

    try:
        models.PDisDataset.objects.get(name=source, version=version, licenced=False).link
    except:
        return True
    return False


def is_licenced_drdi_source(source):
    version = get_today_version()
    try:
        version = get_nedrex_source_version(source)
    except RetryError:
        pass

    try:
        models.DrDiDataset.objects.get(name=source, version=version, licenced=False).link
    except:
        return True
    return False


def remove_old_pdi_data(new_datasets, licenced):
    for dataset in new_datasets:
        try:
            for d in models.PDIDataset.objects.filter(name=dataset.name, licenced=licenced):
                if d != dataset:
                    d.delete()
        except:
            continue


def remove_old_ppi_data(new_datasets, licenced):
    for dataset in new_datasets:
        try:
            for d in models.PPIDataset.objects.filter(name=dataset.name, licenced=licenced):
                if d != dataset:
                    d.delete()
        except:
            continue


def remove_old_pdis_data(new_datasets, licenced):
    for dataset in new_datasets:
        try:
            for d in models.PDisDataset.objects.filter(name=dataset.name, licenced=licenced):
                if d != dataset:
                    d.delete()
        except:
            continue


def remove_old_drdi_data(new_datasets, licenced):
    for dataset in new_datasets:
        try:
            for d in models.DrDiDataset.objects.filter(name=dataset.name, licenced=licenced):
                if d != dataset:
                    d.delete()
        except:
            continue
