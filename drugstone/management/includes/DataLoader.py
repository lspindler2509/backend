import pandas as pd
import json


class DataLoader:
    PATH_PROTEINS = 'data/Proteins/'
    # PATH_DRUGS = 'data/Drugs/'
    PATH_EXPR = 'data/Expression/'
    # PATH_DISORDERS = 'data/Disorders/'
    PATH_PDI = 'data/PDI/'
    PATH_PPI = 'data/PPI/'
    # PATH_PDi = 'data/PDi/'
    PATH_DDi = 'data/DrDi/'

    # Proteins
    # PROTEINS_COVEX = 'protein_list.csv'
    ENTREZ_TO_ENSG = 'entrez_to_ensg.json'

    # Disorders
    # DISORDERS_MONDO = 'disorders.tsv'
    #Drugs
    # DRUG_FILE = 'drug-file.txt'

    #Expressions
    EXPR_FILE = 'gene_tissue_expression.gct'

    # Protein-Protein-Interactions
    PPI_APID = 'apid_9606_Q2.txt'
    # PPI_BIOGRID = 'BIOGRID-ORGANISM-Homo_sapiens-3.5.187.mitab.txt'
    PPI_STRING = 'string_interactions.csv'
    # Protein-Drug-Interactions
    PDI_CHEMBL = 'chembl_drug_gene_interactions_uniq.csv'
    PDI_DGIDB = 'DGIdb_drug_gene_interactions.csv'
    # PDI_DRUGBANK = 'drugbank_drug_gene_interactions_uniq.csv'

    # Protein-Disorder-Interaction
    # PDi_DISGENET = 'disgenet-protein_disorder_association.tsv'

    # Drug-Disorder-Indictation
    DDi_DRUGBANK = 'drugbank-drug_disorder_indication.tsv'

    @staticmethod
    def _clean_entrez(x):
        # convert to string if not empty
        try:
            return str(int(x))
        except Exception:
            return x

    # @staticmethod
    # def _clean_mondo(x):
    #     # convert to string if not empty
    #     try:
    #         return str(int(x))
    #     except Exception:
    #         return x

    # @staticmethod
    # def load_proteins() -> pd.DataFrame:
    #     """Loads the list of proteins used in CoVex
    #
    #     Returns:
    #         pd.DataFrame: columns 'protein_ac', 'gene_name', 'protein_name', 'entrez_id'
    #     """
    #
    #     df = pd.read_csv(f'{DataLoader.PATH_PROTEINS}{DataLoader.PROTEINS_COVEX}')
    #     df['entrez_id'] = df['entrez_id'].map(DataLoader._clean_entrez)
    #     return df

    # @staticmethod
    # def load_drugs()-> pd.DataFrame:
    #     return pd.read_csv(f'{DataLoader.PATH_DRUGS}{DataLoader.DRUG_FILE}', sep='\t')

    @staticmethod
    def load_expressions() -> pd.DataFrame:
        return pd.read_csv(f'{DataLoader.PATH_EXPR}{DataLoader.EXPR_FILE}', sep='\t')


    # @staticmethod
    # def load_disorders() -> pd.DataFrame:
    #     """Loads the list of disorders used in Nedrex
    #
    #     Returns:
    #         pd.DataFrame: columns 'mondo_id', 'disorder_name', 'icd_10'
    #     """
    #
    #     df = pd.read_csv(f'{DataLoader.PATH_DISORDERS}{DataLoader.DISORDERS_MONDO}', sep='\t')
    #     df['mondo_id'] = df['mondo_id'].map(DataLoader._clean_mondo)
    #     return df

    @staticmethod
    def load_ensg() -> dict:
        """Loads the json of entrez -> list of ensg relations for all proteins used in CoVex

        Returns:
            dict with {entrez1: [ensg1, ensg2], ...}
        """

        f = open(f'{DataLoader.PATH_PROTEINS}{DataLoader.ENTREZ_TO_ENSG}', )
        data = json.load(f)
        return data

    @staticmethod
    def load_ppi_string() -> pd.DataFrame:
        """Loads the STRING PPI interactions with Entrez IDs

        Returns:
            pd.DataFrame: columns 'entrez_a', 'entrez_b'
        """
        df = pd.read_csv(f'{DataLoader.PATH_PPI}{DataLoader.PPI_STRING}', index_col=0)
        df['entrez_a'] = df['entrez_a'].map(DataLoader._clean_entrez)
        df['entrez_b'] = df['entrez_b'].map(DataLoader._clean_entrez)
        df = df.drop_duplicates()
        return df

    # @staticmethod
    # def load_ppi_biogrid() -> pd.DataFrame:
    #     """Loads the Biogrid PPI interactions with Entex IDs
    #
    #     Returns:
    #         pd.DataFrame: columns 'entrez_a', 'entrez_b'
    #     """
    #     df = pd.read_csv(f'{DataLoader.PATH_PPI}{DataLoader.PPI_BIOGRID}', sep='\t')[
    #         ['#ID Interactor A', 'ID Interactor B', 'Interaction Detection Method', 'Publication Identifiers',
    #          'Interaction Types', 'Confidence Values']]
    #
    #     def parse_interactor(x):
    #         # format entrez protein/locuslink:6416
    #         # wanted: 6416
    #         return x.split(':')[-1]
    #
    #     def parse_publication_identifiers(x):
    #         # format: pubmed:9006895
    #         # wanted: 9006895
    #         return x.split(':')[-1]
    #
    #     def parse_interaction_detection_method_get_id(x):
    #         # format: psi-mi:"MI:0018"(two hybrid)
    #         # wanted: MI:0018
    #         return x.split('"')[1]
    #
    #     def parse_interaction_detection_method_get_name(x):
    #         # format: psi-mi:"MI:0018"(two hybrid)
    #         # wanted: two hybrid
    #         return x.split('(')[1][:-1]
    #
    #     def parse_interaction_types_get_id(x):
    #         # format: psi-mi:"MI:0407"(direct interaction)
    #         # wanted: MI:0407
    #         return x.split('"')[1]
    #
    #     def parse_interaction_types_get_name(x):
    #         # format: psi-mi:"MI:0407"(direct interaction)
    #         # wanted: direct interaction
    #         return x.split('(')[1][:-1]
    #
    #     def parse_confidence_value(x):
    #         # format: score:7.732982515 or '-'
    #         # wanted: 7.732982515 or '-'
    #         if x == '-':
    #             return '-'
    #         else:
    #             return x.split(':')[1]
    #
    #     df['#ID Interactor A'] = df['#ID Interactor A'].map(parse_interactor)
    #     df['ID Interactor B'] = df['ID Interactor B'].map(parse_interactor)
    #
    #     df['Interaction Detection Method ID'] = df['Interaction Detection Method'].map(
    #         parse_interaction_detection_method_get_id)
    #     df['Interaction Detection Method Name'] = df['Interaction Detection Method'].map(
    #         parse_interaction_detection_method_get_name)
    #     df = df.drop('Interaction Detection Method', axis=1)
    #
    #     df['Publication Identifiers'] = df['Publication Identifiers'].map(parse_publication_identifiers)
    #
    #     df['Interaction Types ID'] = df['Interaction Types'].map(parse_interaction_types_get_id)
    #     df['Interaction Types Name'] = df['Interaction Types'].map(parse_interaction_types_get_name)
    #     df = df.drop('Interaction Types', axis=1)
    #
    #     df['Confidence Values'] = df['Confidence Values'].map(parse_confidence_value)
    #
    #     # remove dirty data (entrez id is sometimes protein name
    #     to_remove = ['P0DTC1', 'P0DTD2', 'Q7TLC7']
    #     df = df[~df['#ID Interactor A'].isin(to_remove)]
    #
    #     df = df.rename(
    #         columns={
    #             '#ID Interactor A': 'entrez_a',
    #             'ID Interactor B': 'entrez_b',
    #             'Publication Identifiers': 'pubmed_id',
    #             'Confidence Values': 'confidence_value',
    #             'Interaction Detection Method ID': 'detection_method_psi_mi',
    #             'Interaction Detection Method Name': 'detection_method_name',
    #             'Interaction Types ID': 'type_psi_mi',
    #             'Interaction Types Name': 'type_name'
    #         }
    #     )
    #
    #     df['entrez_a'] = df['entrez_a'].map(DataLoader._clean_entrez)
    #     df['entrez_b'] = df['entrez_b'].map(DataLoader._clean_entrez)
    #
    #     return df

    @staticmethod
    def load_ppi_apid() -> pd.DataFrame:
        """Loads the APID PPI interactions with Uniprot ACs

        Returns:
            pd.DataFrame: columns 'from_protein_ac', 'to_protein_ac'
        """
        df = pd.read_csv(f'{DataLoader.PATH_PPI}{DataLoader.PPI_APID}', index_col=0, sep='\t')
        df = df.rename(columns={
            'UniprotID_A': 'from_protein_ac',
            'UniprotID_B': 'to_protein_ac'
        })
        return df[['from_protein_ac', 'to_protein_ac']]

    @staticmethod
    def load_pdi_chembl() -> pd.DataFrame:
        """Loads the Chembl PDI interactions with DrugBank drug IDs and Uniprot ACs

        Returns:
            pd.DataFrame: columns "drug_id" and "protein_ac"
        """
        return pd.read_csv(f'{DataLoader.PATH_PDI}{DataLoader.PDI_CHEMBL}')

    # @staticmethod
    # def load_pdis_disgenet() -> pd.DataFrame:
    #     """Loads the DisGeNET PDis associations with UniprotAC Numbers and Mondo IDs
    #
    #     Returns:
    #         pd.DataFrame: columns "protein_name", "disorder_name" and "score"
    #     """
    #     return pd.read_csv(f'{DataLoader.PATH_PDi}{DataLoader.PDi_DISGENET}', sep='\t', dtype={'disorder_name':str, 'protein_name':str, 'score':float})

    @staticmethod
    def load_drdis_drugbank() -> pd.DataFrame:
        """Loads the DrugBank DrDi indications with DrugBank and Mondo IDs

        Returns:
            pd.DataFrame: columns "drugbank_id" and "mondo_id"
        """
        return pd.read_csv(f'{DataLoader.PATH_DDi}{DataLoader.DDi_DRUGBANK}', sep='\t', dtype={'drugbank_id':str, 'mondo_id':str})

    @staticmethod
    def load_pdi_dgidb() -> pd.DataFrame:
        """Loads the DGIdb PDI interactions with DrugBank drug IDs and Entrez Gene IDs

        Returns:
            pd.DataFrame: columns "drug_id" and "entrez_id"
        """
        df = pd.read_csv(f'{DataLoader.PATH_PDI}{DataLoader.PDI_DGIDB}', index_col=0)
        df['entrez_id'] = df['entrez_id'].map(DataLoader._clean_entrez)
        return df

    # @staticmethod
    # def load_pdi_drugbank() -> pd.DataFrame:
    #     """Loads the drugbank PDI interactions with DrugBank drug IDs and Entrez Gene IDs
    #
    #     Returns:
    #         pd.DataFrame: columns "drug_id" and "entrez_id"
    #     """
    #     df = pd.read_csv(f'{DataLoader.PATH_PDI}{DataLoader.PDI_DRUGBANK}').dropna()
    #     df['entrez_id'] = df['entrez_id'].map(DataLoader._clean_entrez)
    #     return df
