from django.core.management.base import BaseCommand
from django.db import OperationalError

from drugstone.models import Protein, Drug, Tissue, ExpressionLevel, PPIDataset, PDIDataset, Disorder, PDisDataset, \
    DrDiDataset, EnsemblGene, CellularComponent
from drugstone.models import ProteinProteinInteraction, ProteinDrugInteraction, ProteinDisorderAssociation, \
    DrugDisorderIndication, ActiveIn

from drugstone.management.includes.DataPopulator import DataPopulator
from .import_from_nedrex import NedrexImporter
from drugstone.management.includes.NodeCache import NodeCache
from drugstone.management.includes import DatasetLoader
# from ..includes.DatasetLoader import remove_old_pdi_data, remove_old_ppi_data, remove_old_pdis_data, \
#     remove_old_drdi_data


class DatabasePopulator:
    def __init__(self, data_dir):
        self.data_dir = data_dir

    def delete_model(self, model):
        from django.db import connection
        cursor = connection.cursor()
        try:
            cursor.execute('TRUNCATE TABLE "{0}" CASCADE'.format(model._meta.db_table))
        except OperationalError:
            cursor.execute('DELETE FROM "{0}"'.format(model._meta.db_table))

    def delete_all(self):
        models = ['PPI', 'PDI', 'DrDi', 'Protein', 'Drug', 'Disorder', 'PDi', 'Expression', 'Tissue', 'CellularComponent']
        self.delete_models(models)

    def delete_models(self, model_list):
        for model_name in model_list:
            print(f'Deleting {model_name} model ...')

            if model_name == 'PPI':
                self.delete_model(PPIDataset)
                self.delete_model(ProteinProteinInteraction)
            elif model_name == 'PDI':
                self.delete_model(PDIDataset)
                self.delete_model(ProteinDrugInteraction)
            elif model_name == 'DrDi':
                self.delete_model(DrDiDataset)
                self.delete_model(DrugDisorderIndication)
            elif model_name == 'Protein':
                self.delete_model(Protein)
                self.delete_model(EnsemblGene)
            elif model_name == 'Drug':
                self.delete_model(Drug)
            elif model_name == 'Disorder':
                self.delete_model(Disorder)
            elif model_name == 'PDi':
                self.delete_model(PDisDataset)
                self.delete_model(ProteinDisorderAssociation)
            elif model_name == 'Expression':
                self.delete_model(ExpressionLevel)
            elif model_name == 'Tissue':
                self.delete_model(Tissue)
            elif model_name == 'CellularComponent':
                self.delete_model(CellularComponent)
                self.delete_model(ActiveIn)


class Command(BaseCommand):
    def add_arguments(self, parser):
        # dataset directory
        parser.add_argument('-dd', '--data_dir', type=str, help='Dataset directory path')
        parser.add_argument('-dm', '--delete_model', type=str, help='Delete model(s)')
        parser.add_argument('-c', '--clear', action='store_true', help='Delete all models')
        parser.add_argument('-a', '--all', action='store_true', help='Populate all tables')
        parser.add_argument('-u', '--update', action='store_true', help='Execute database update for selected tables')

        parser.add_argument('-p', '--proteins', action='store_true', help='Populate Proteins')
        parser.add_argument('-di', '--disorders', action='store_true', help='Populate Disorders')
        parser.add_argument('-dr', '--drugs', action='store_true', help='Drug file name')
        
        parser.add_argument('-cc', '--cellular_components', action='store_true', help='Populate cellular components')

        parser.add_argument('-exp', '--exp', action='store_true',
                            help='Tissue expression file (.gct without first 2 lines)')

        parser.add_argument('-pp', '--protein_protein', action='store_true',
                            help='Populate Protein-Protein Interactions')
        parser.add_argument('-pdr', '--protein_drug', action='store_true', help='Populate Protein-Drug Interactions')
        parser.add_argument('-pdi', '--protein_disorder', action='store_true',
                            help='Populate Protein-Disorder Associations')
        parser.add_argument('-ddi', '--drug_disorder', action='store_true', help='Populate Drug-Disorder Indications')
        parser.add_argument('-t', '--test', action='store_true', help='Running some function on startup')

    def handle(self, *args, **kwargs):
        populate(kwargs)


def populate(kwargs):
    nedrex_api_url_open = "https://api.nedrex.net/open"
    nedrex_api_url_licensed = "https://api.nedrex.net/licensed"

    data_dir = kwargs['data_dir']

    db_populator = DatabasePopulator(data_dir=data_dir)
    if 'test' in kwargs and kwargs['test']:
        pass
        # remove_old_ppi_data([PPIDataset.objects.filter(name='biogrid', licenced=False).last()], False)
        # remove_old_ppi_data([PPIDataset.objects.filter(name='iid', licenced=False).last()], False)
        # remove_old_ppi_data([PPIDataset.objects.filter(name='intact', licenced=False).last()], False)
        # remove_old_pdis_data([PDisDataset.objects.filter(name='disgenet', licenced=False).last()], False)
        # remove_old_pdis_data([PDisDataset.objects.filter(name='omim', licenced=True).last()], True)
        # remove_old_drdi_data([DrDiDataset.objects.filter(name='ctd', licenced=False).last()], False)
        # remove_old_drdi_data([DrDiDataset.objects.filter(name='drugcentral', licenced=False).last()], False)
    if 'clear' in kwargs and kwargs['clear']:
        db_populator.delete_all()

    if 'delete_model' in kwargs and kwargs['delete_model'] is not None:
        model_list = kwargs['delete_model'].split(',')
        db_populator.delete_models(model_list)

    cache = NodeCache()
    update = True if kwargs['update'] else False
    importer = NedrexImporter(nedrex_api_url_licensed, nedrex_api_url_open, cache)
    populator = DataPopulator(cache)

    total_n = 0
    nedrex_update = False


    if 'all' in kwargs and kwargs['all']:
        kwargs['drugs'] = True
        kwargs['disorders'] = True
        kwargs['proteins'] = True
        kwargs['exp'] = True
        kwargs['protein_protein'] = True
        kwargs['protein_drug'] = True
        kwargs['protein_disorder'] = True
        kwargs['drug_disorder'] = True
        kwargs['cellular_components'] = True

    if kwargs['drugs']:
        print('Populating Drugs...')
        n = NedrexImporter.import_drugs(importer, update)
        total_n += n
        nedrex_update = True
        print(f'Populated {n} Drugs.')

    if kwargs['disorders']:
        print('Populating Disorders...')
        n = NedrexImporter.import_disorders(importer, update)
        total_n += n
        nedrex_update = True
        print(f'Populated {n} Disorders.')

    if kwargs['proteins']:
        print('Populating Proteins...')
        n = NedrexImporter.import_proteins(importer, update)
        total_n += n
        nedrex_update = True
        print(f'Populated {n} Proteins.')
        print('Populating ENSG IDs...')
        n = DataPopulator.populate_ensg(populator, update)
        total_n += n
        print(f'Populated {n} ENSG IDs.')

    if kwargs['exp']:
        print('Populating Expressions...')
        n = DataPopulator.populate_expressions(populator, update)
        total_n += n
        print(f'Populated {n} Expressions.')

    if kwargs['protein_drug']:
        print('Importing PDIs from unlicensed NeDRexDB...')
        n = NedrexImporter.import_drug_target_interactions(importer,
                                                           DatasetLoader.get_drug_target_nedrex(nedrex_api_url_open,
                                                                                                False),
                                                           update)
        total_n += n
        print(f'Imported {n} PDIs from unlicensed NeDRexDB')

        print('Importing PDIs from licensed NeDRexDB...')
        n = NedrexImporter.import_drug_target_interactions(importer,
                                                           DatasetLoader.get_drug_target_nedrex(nedrex_api_url_licensed,
                                                                                                True),
                                                           update)
        total_n += n
        nedrex_update = True
        print(f'Imported {n} PDIs from licensed NeDRexDB')

        dataset, created = DatasetLoader.get_drug_target_chembl()
        if created:
            print('Populating PDIs from Chembl...')
            n = DataPopulator.populate_pdi_chembl(populator, dataset, update)
            total_n += n
            print(f'Populated {n} PDIs from Chembl.')
        else:
            print('Chembl already populated.')

        dataset, created = DatasetLoader.get_drug_target_dgidb()
        if created:
            print('Populating PDIs from DGIdb...')
            n = DataPopulator.populate_pdi_dgidb(populator, dataset, update)
            total_n += n
            print(f'Populated {n} PDIs from DGIdb.')
        else:
            print('DGIdb already populated.')

    if kwargs['protein_disorder']:
        print('Importing PDis from unlicensed NeDRexDB...')
        n = NedrexImporter.import_protein_disorder_associations(importer,
                                                                DatasetLoader.get_protein_disorder_nedrex(
                                                                    nedrex_api_url_open, False),
                                                                update)
        total_n += n
        print(f'Imported {n} PDis from unlicensed NeDRexDB')

        print('Importing PDis from licenced NeDRexDB...')
        n = NedrexImporter.import_protein_disorder_associations(importer,
                                                                DatasetLoader.get_protein_disorder_nedrex(
                                                                    nedrex_api_url_licensed, True),
                                                                update)
        total_n += n
        nedrex_update = True
        print(f'Imported {n} PDis from licenced NeDRexDB')

    if kwargs['drug_disorder']:
        print('Importing DrDis from unlicensed NeDRexDB...')
        n = NedrexImporter.import_drug_disorder_indications(importer,
                                                            DatasetLoader.get_drug_disorder_nedrex(nedrex_api_url_open,
                                                                                                   False),
                                                            update)
        total_n += n
        print(f'Imported {n} DrDis from unlicensed NeDRexDB')

        print('Importing DrDis from licenced NeDRexDB...')
        n = NedrexImporter.import_drug_disorder_indications(importer,
                                                            DatasetLoader.get_drug_disorder_nedrex(
                                                                nedrex_api_url_licensed, True),
                                                            update)
        total_n += n
        nedrex_update = True
        print(f'Imported {n} DrDis from licenced NeDRexDB')
        print('Populating DrDi indications from DrugBank...')
        n = DataPopulator.populate_drdis_drugbank(populator, DatasetLoader.get_drug_disorder_drugbank(), update)
        total_n += n
        print(f'Populated {n} DrDi associations from DrugBank.')

    if kwargs['cellular_components']:
        print('Importing cellular components...')
        n = NedrexImporter.import_cellularComponent(populator, update)
        print(f'Populated {n} Cellular components relations.')
        
    
    if kwargs['protein_protein']:
        print('Importing PPIs from unlicensed NeDRexDB...')
        n = NedrexImporter.import_protein_protein_interactions(importer,
                                                               DatasetLoader.get_ppi_nedrex(nedrex_api_url_open, False),
                                                               update)
        total_n += n
        print(f'Imported {n} PPIs from unlicensed NeDRexDB')
        print('Importing PPIs from licenced NeDRexDB...')
        n = NedrexImporter.import_protein_protein_interactions(importer,
                                                               DatasetLoader.get_ppi_nedrex(nedrex_api_url_licensed,
                                                                                            True),
                                                               update)
        total_n += n
        nedrex_update = True
        print(f'Imported {n} PPIs from licensed NeDRexDB')
        
        dataset, created = DatasetLoader.get_ppi_string()
        if created:
            print('Populating PPIs from STRING...')
            n = DataPopulator.populate_ppi_string(populator, dataset, update)
            total_n += n
            print(f'Populated {n} PPIs from STRING.')
        else:
            print('STRING already populated.')

        dataset, created = DatasetLoader.get_ppi_apid()
        if created:
            print('Populating PPIs from APID...')
            n = DataPopulator.populate_ppi_apid(populator, dataset, update)
            total_n += n
            print(f'Populated {n} PPIs from APID.')
        else:
            print('APID already populated.')

    if nedrex_update:
        from drugstone.management.includes.DatasetLoader import update_license
        update_license()

    cache.clear()
    return total_n
