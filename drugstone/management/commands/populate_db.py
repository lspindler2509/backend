from django.core.management.base import BaseCommand
from django.db import OperationalError

from drugstone.models import Protein, Drug, Tissue, ExpressionLevel, PPIDataset, PDIDataset, Disorder, PDisDataset, \
    DrDiDataset, EnsemblGene
from drugstone.models import ProteinProteinInteraction, ProteinDrugInteraction, ProteinDisorderAssociation, \
    DrugDisorderIndication

from drugstone.management.includes.DataPopulator import DataPopulator
from .import_from_nedrex import NedrexImporter
from drugstone.management.includes.NodeCache import NodeCache
from drugstone.management.includes import DatasetLoader




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
        models = ['PPI', 'PDI', 'DrDi', 'Protein', 'Drug', 'Disorder', 'PDi', 'Expression', 'Tissue']
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

        parser.add_argument('-exp', '--exp', action='store_true',
                            help='Tissue expression file (.gct without first 2 lines)')

        parser.add_argument('-pp', '--protein_protein', action='store_true',
                            help='Populate Protein-Protein Interactions')
        parser.add_argument('-pdr', '--protein_drug', action='store_true', help='Populate Protein-Drug Interactions')
        parser.add_argument('-pdi', '--protein_disorder', action='store_true',
                            help='Populate Protein-Disorder Associations')
        parser.add_argument('-ddi', '--drug_disorder', action='store_true', help='Populate Drug-Disorder Indications')

    def handle(self, *args, **kwargs):
        populate(kwargs)

def populate(kwargs):

    nedrex_api_url = "http://82.148.225.92:8123/"
    data_dir = kwargs['data_dir']

    db_populator = DatabasePopulator(data_dir=data_dir)

    if kwargs['clear']:
        db_populator.delete_all()

    if kwargs['delete_model'] is not None:
        model_list = kwargs['delete_model'].split(',')
        db_populator.delete_models(model_list)

    cache = NodeCache()
    update = True if kwargs['update'] else False
    importer = NedrexImporter(nedrex_api_url, cache)
    populator = DataPopulator(cache)

    if kwargs['all']:
        kwargs['drugs'] = True
        kwargs['disorders'] = True
        kwargs['proteins'] = True
        kwargs['exp'] = True
        kwargs['protein_protein'] = True
        kwargs['protein_drug'] = True
        kwargs['protein_disorder'] = True
        kwargs['drug_disorder'] = True

    if kwargs['drugs']:
        print('Populating Drugs...')
        n = NedrexImporter.import_drugs(importer, update)
        print(f'Populated {n} Drugs.')

    if kwargs['disorders']:
        print('Populating Disorders...')
        n = NedrexImporter.import_disorders(importer, update)
        print(f'Populated {n} Disorders.')

    if kwargs['proteins']:
        print('Populating Proteins...')
        n = NedrexImporter.import_proteins(importer, update)
        print(f'Populated {n} Proteins.')
        print('Populating ENSG IDs...')
        n = DataPopulator.populate_ensg(populator, update)
        print(f'Populated {n} ENSG IDs.')

    if kwargs['exp']:
        print('Populating Expressions...')
        n = DataPopulator.populate_expressions(populator, update)
        print(f'Populated {n} Expressions.')

    if kwargs['protein_protein']:
        print('Importing PPIs from NeDRexDB...')
        n = NedrexImporter.import_protein_protein_interactions(importer,
                                                               DatasetLoader.get_ppi_nedrex(nedrex_api_url),
                                                               update)
        print(f'Imported {n} PPIs from NeDRexDB')
        print('Populating PPIs from STRING...')
        n = DataPopulator.populate_ppi_string(populator, DatasetLoader.get_ppi_string(), update)
        print(f'Populated {n} PPIs from STRING.')

        print('Populating PPIs from APID...')
        n = DataPopulator.populate_ppi_apid(populator, DatasetLoader.get_ppi_apid(), update)
        print(f'Populated {n} PPIs from APID.')

        print('Populating PPIs from BioGRID...')
        n = DataPopulator.populate_ppi_biogrid(populator, DatasetLoader.get_ppi_biogrid(), update)
        print(f'Populated {n} PPIs from BioGRID.')

    if kwargs['protein_drug']:
        print('Importing PDIs from NeDRexDB...')
        n = NedrexImporter.import_drug_target_interactions(importer,
                                                           DatasetLoader.get_drug_target_nedrex(nedrex_api_url),
                                                           update)
        print(f'Imported {n} PDIs from NeDRexDB')

        print('Populating PDIs from Chembl...')
        n = DataPopulator.populate_pdi_chembl(populator, DatasetLoader.get_drug_target_chembl(), update)
        print(f'Populated {n} PDIs from Chembl.')

        print('Populating PDIs from DGIdb...')
        n = DataPopulator.populate_pdi_dgidb(populator, DatasetLoader.get_drug_target_dgidb(), update)
        print(f'Populated {n} PDIs from DGIdb.')

        print('Populating PDIs from DrugBank...')
        n = DataPopulator.populate_pdi_drugbank(populator, DatasetLoader.get_drug_target_drugbank(), update)
        print(f'Populated {n} PDIs from DrugBank.')

    if kwargs['protein_disorder']:
        print('Importing PDis from NeDRexDB...')
        n = NedrexImporter.import_protein_disorder_associations(importer,
                                                                DatasetLoader.get_protein_disorder_nedrex(
                                                                    nedrex_api_url),
                                                                update)
        print(f'Imported {n} PDis from NeDRexDB')
        print('Populating PDis associations from DisGeNET...')
        n = DataPopulator.populate_pdis_disgenet(populator, DatasetLoader.get_disorder_protein_disgenet(), update)
        print(f'Populated {n} PDis associations from DisGeNET.')

    if kwargs['drug_disorder']:
        print('Importing DrDis from NeDRexDB...')
        n = NedrexImporter.import_drug_disorder_indications(importer,
                                                            DatasetLoader.get_drug_disorder_nedrex(nedrex_api_url),
                                                            update)
        print(f'Imported {n} DrDis from NeDRexDB')
        print('Populating DrDi indications from DrugBank...')
        n = DataPopulator.populate_drdis_drugbank(populator, DatasetLoader.get_drug_disorder_drugbank(), update)
        print(f'Populated {n} DrDi associations from DrugBank.')
