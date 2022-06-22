from django.core.management.base import BaseCommand
import pandas as pd
from django.db import OperationalError, IntegrityError

from drugstone.models import Protein, Drug, Tissue, ExpressionLevel, PPIDataset, PDIDataset, Disorder, PDisDataset
from drugstone.models import ProteinProteinInteraction, ProteinDrugInteraction

from drugstone.management.includes.DataPopulator import DataPopulator


class DatabasePopulator:
    def __init__(self, data_dir,
                #  protein_file,
                 drug_file,
                #  protein_protein_interaction_file,
                #  protein_drug_interaction_file,
                 tissue_expression_file,
                 ):
        self.data_dir = data_dir
        # self.protein_file = protein_file
        self.drug_file = drug_file
        # self.ppi_file = protein_protein_interaction_file
        # self.pdi_file = protein_drug_interaction_file
        self.exp_file = tissue_expression_file

    def delete_model(self, model):
        from django.db import connection
        cursor = connection.cursor()
        try:
            cursor.execute('TRUNCATE TABLE "{0}" CASCADE'.format(model._meta.db_table))
        except OperationalError:
            cursor.execute('DELETE FROM "{0}"'.format(model._meta.db_table))

    def delete_models(self, model_list):
        for model_name in model_list:
            print(f'Deleting {model_name} model ...')

            if model_name == 'PPI':
                self.delete_model(ProteinProteinInteraction)
            elif model_name == 'PDI':
                self.delete_model(ProteinDrugInteraction)
            elif model_name == 'Protein':
                self.delete_model(Protein)
            elif model_name == 'Drug':
                self.delete_model(Drug)
            elif model_name == 'Disorder':
                self.delete_model(Disorder)
            elif model_name == 'PDiAssociations':
                self.delete_model(PDisDataset)
            elif model_name == 'Tissue':
                self.delete_model(Tissue)
            elif model_name == 'PPIDataset':
                self.delete_model(PPIDataset)
            elif model_name == 'PDIDataset':
                self.delete_model(PDIDataset)
            elif model_name == 'DrDiDataset':
                self.delete_model(DrDiDataset)


    def populate_drug_model(self):
        print('Populating Drug model ...')
        drug_df = pd.read_csv(f'{self.data_dir}/data_drugstone/{self.drug_file}', delimiter='\t')
        for _, row in drug_df.iterrows():
            drug_id = row['drug_id']
            drug_name = row['drug_name']
            drug_status = row['drug_status']
            # links = row['links']
            Drug.objects.update_or_create(
                drug_id=drug_id, 
                name=drug_name, 
                status=drug_status, 
                # links=links
                )

        print('Done!\n')


    def populate_exp_model(self):
        print('Populating Tissue and ExpressionLevel model ...')
        exp_df = pd.read_csv(f'{self.data_dir}/data_drugstone/{self.exp_file}', delimiter='\t')

        tissues_models = dict()
        for tissue_name in exp_df.columns.values[2:]:
            try:
                tissue_model = Tissue.objects.get(name=tissue_name)
            except Tissue.DoesNotExist:
                tissue_model = Tissue.objects.create(name=tissue_name)
            tissues_models[tissue_name] = tissue_model

        proteins_linked = 0

        for _, row in exp_df.iterrows():
            gene_name = row['Description']

            for protein_model in Protein.objects.filter(gene=gene_name).all():
                proteins_linked += 1

                for tissue_name, tissue_model in tissues_models.items():
                    try:
                        ExpressionLevel.objects.create(protein=protein_model,
                                                       tissue=tissue_model,
                                                       expression_level=row[tissue_name])
                    except IntegrityError:
                        pass

        print(f'Added {proteins_linked} expression levels!\n')



class Command(BaseCommand):
    def add_arguments(self, parser):

        # dataset directory
        parser.add_argument('-dd', '--data_dir', type=str, help='Dataset directory path')
        # parser.add_argument('-p', '--protein_file', type=str, help='Protein file')
        parser.add_argument('-dr', '--drug_file', type=str, help='Drug file name')
        # parser.add_argument('-ppi', '--ppi_file', type=str, help='Protein-Protein interaction file')
        # parser.add_argument('-pdi', '--pdi_file', type=str, help='Protein-Drug interaction file')
        parser.add_argument('-exp', '--exp_file', type=str, help='Tissue expression file (.gct without first 2 lines)')
        parser.add_argument('-dm', '--delete_model', type=str, help='Delete model(s)')

        parser.add_argument('-p', '--proteins', type=str, help='Populate Proteins')
        parser.add_argument('-di', '--disorders', type=str, help='Populate Disorders')
        parser.add_argument('-pp', '--protein_protein', type=str, help='Populate Protein-Protein Interactions')
        parser.add_argument('-pdr', '--protein_drug', type=str, help='Populate Protein-Drug Interactions')
        parser.add_argument('-pdi', '--protein_disorder', type=str, help='Populate Protein-Disorder Associations')
        parser.add_argument('-ddi', '--drug_disorder', type=str, help='Populate Drug-Disorder Indications')

    def handle(self, *args, **kwargs):

        data_dir = kwargs['data_dir']
        # protein_file = kwargs['protein_file']
        drug_file = kwargs['drug_file']
        # ppi_file = kwargs['ppi_file']
        # pdi_file = kwargs['pdi_file']
        exp_file = kwargs['exp_file']

        p = kwargs['proteins']
        pp = kwargs['protein_protein']
        pd = kwargs['protein_drug']

        db_populator = DatabasePopulator(data_dir=data_dir,
                                        # protein_file=protein_file,
                                        drug_file=drug_file,
                                        # protein_protein_interaction_file=ppi_file,
                                        # protein_drug_interaction_file=pdi_file,
                                        tissue_expression_file=exp_file,
                                        )

        if kwargs['delete_model'] is not None:
            model_list = kwargs['delete_model'].split(',')
            db_populator.delete_models(model_list)
            return

        if kwargs['drug_file'] is not None:
            db_populator.populate_drug_model()

        # if kwargs['protein_file'] is not None:
        #     db_poulator.populate_protein_model()

        # if kwargs['pdi_file'] is not None:
        #     db_poulator.populate_pdi_model()

        # if kwargs['ppi_file'] is not None:
        #     db_poulator.populate_ppi_model()

        if kwargs['exp_file'] is not None:
            db_populator.populate_exp_model()

        if kwargs['proteins'] is not None:
            print('Populating Proteins...')
            n = DataPopulator.populate_proteins()
            print(f'Populated {n} Proteins.')
            
            print('Populating ENSG IDs...')
            n = DataPopulator.populate_ensg()
            print(f'Populated {n} ENSG IDs.')

        if kwargs['disorders'] is not None:
            print('Populating Disorders...')
            n = DataPopulator.populate_disorders()
            print(f'Populated {n} Disorders.')

        if kwargs['protein_protein'] is not None:
            print('Populating PPIs from STRING...')
            n = DataPopulator.populate_ppi_string()
            print(f'Populated {n} PPIs from STRING.')

            print('Populating PPIs from APID...')
            n = DataPopulator.populate_ppi_apid()
            print(f'Populated {n} PPIs from APID.')

            print('Populating PPIs from BioGRID...')
            n = DataPopulator.populate_ppi_biogrid()
            print(f'Populated {n} PPIs from BioGRID.')

        if kwargs['protein_drug'] is not None:
            print('Populating PDIs from Chembl...')
            n = DataPopulator.populate_pdi_chembl()
            print(f'Populated {n} PDIs from Chembl.')

            print('Populating PDIs from DGIdb...')
            n = DataPopulator.populate_pdi_dgidb() 
            print(f'Populated {n} PDIs from DGIdb.')

            print('Populating PDIs from DrugBank...')
            n = DataPopulator.populate_pdi_drugbank()
            print(f'Populated {n} PDIs from DrugBank.')
        if kwargs['protein_disorder'] is not None:
            print('Populating PDis associations from DisGeNET...')
            n=DataPopulator.populate_pdis_disgenet()
            print(f'Populated {n} PDis associations from DisGeNET.')
        if kwargs['drug_disorder'] is not None:
            print('Populating DrDi indications from DrugBank...')
            n=DataPopulator.populate_drdis_drugbank()
            print(f'Populated {n} DrDi associations from DrugBank.')
