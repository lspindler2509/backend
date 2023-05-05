import pandas as pd
import json


class DataLoader:
    PATH_PROTEINS = "data/Proteins/"
    # PATH_DRUGS = 'data/Drugs/'
    PATH_EXPR = "data/Expression/"
    # PATH_DISORDERS = 'data/Disorders/'
    PATH_PDI = "data/PDI/"
    PATH_PPI = "data/PPI/"
    # PATH_PDi = 'data/PDi/'
    PATH_DDi = "data/DrDi/"

    # Proteins
    # PROTEINS_COVEX = 'protein_list.csv'
    ENTREZ_TO_ENSG = "entrez_to_ensg.json"

    # Disorders
    # DISORDERS_MONDO = 'disorders.tsv'
    # Drugs
    # DRUG_FILE = 'drug-file.txt'

    # Expressions
    EXPR_FILE = "gene_tissue_expression.gct"

    # Protein-Protein-Interactions
    PPI_APID = "apid_9606_Q2.txt"
    # PPI_BIOGRID = 'BIOGRID-ORGANISM-Homo_sapiens-3.5.187.mitab.txt'
    PPI_STRING = "string_interactions.csv"
    # Protein-Drug-Interactions
    PDI_CHEMBL = "chembl_drug_gene_interactions_uniq.csv"
    PDI_DGIDB = "DGIdb_drug_gene_interactions.csv"
    # PDI_DRUGBANK = 'drugbank_drug_gene_interactions_uniq.csv'

    # Protein-Disorder-Interaction
    # PDi_DISGENET = 'disgenet-protein_disorder_association.tsv'

    # Drug-Disorder-Indictation
    DDi_DRUGBANK = "drugbank-drug_disorder_indication.tsv"

    @staticmethod
    def _clean_entrez(x):
        # convert to string if not empty
        try:
            return str(int(x))
        except Exception:
            return x

    @staticmethod
    def load_expressions() -> pd.DataFrame:
        return pd.read_csv(f"{DataLoader.PATH_EXPR}{DataLoader.EXPR_FILE}", sep="\t")

    @staticmethod
    def load_ensg() -> dict:
        """Loads the json of entrez -> list of ensg relations for all proteins used in CoVex

        Returns:
            dict with {entrez1: [ensg1, ensg2], ...}
        """

        f = open(
            f"{DataLoader.PATH_PROTEINS}{DataLoader.ENTREZ_TO_ENSG}",
        )
        data = json.load(f)
        return data

    @staticmethod
    def load_ppi_string() -> pd.DataFrame:
        """Loads the STRING PPI interactions with Entrez IDs

        Returns:
            pd.DataFrame: columns 'entrez_a', 'entrez_b'
        """
        df = pd.read_csv(f"{DataLoader.PATH_PPI}{DataLoader.PPI_STRING}", index_col=0)
        df["entrez_a"] = df["entrez_a"].map(DataLoader._clean_entrez)
        df["entrez_b"] = df["entrez_b"].map(DataLoader._clean_entrez)
        df = df.drop_duplicates()
        return df

    @staticmethod
    def load_ppi_apid() -> pd.DataFrame:
        """Loads the APID PPI interactions with Uniprot ACs

        Returns:
            pd.DataFrame: columns 'from_protein_ac', 'to_protein_ac'
        """
        df = pd.read_csv(
            f"{DataLoader.PATH_PPI}{DataLoader.PPI_APID}", index_col=0, sep="\t"
        )
        df = df.rename(
            columns={"UniprotID_A": "from_protein_ac", "UniprotID_B": "to_protein_ac"}
        )
        return df[["from_protein_ac", "to_protein_ac"]]

    @staticmethod
    def load_pdi_chembl() -> pd.DataFrame:
        """Loads the Chembl PDI interactions with DrugBank drug IDs and Uniprot ACs

        Returns:
            pd.DataFrame: columns "drug_id" and "protein_ac"
        """
        return pd.read_csv(f"{DataLoader.PATH_PDI}{DataLoader.PDI_CHEMBL}")

    @staticmethod
    def load_drdis_drugbank() -> pd.DataFrame:
        """Loads the DrugBank DrDi indications with DrugBank and Mondo IDs

        Returns:
            pd.DataFrame: columns "drugbank_id" and "mondo_id"
        """
        return pd.read_csv(
            f"{DataLoader.PATH_DDi}{DataLoader.DDi_DRUGBANK}",
            sep="\t",
            dtype={"drugbank_id": str, "mondo_id": str},
        )

    @staticmethod
    def load_pdi_dgidb() -> pd.DataFrame:
        """Loads the DGIdb PDI interactions with DrugBank drug IDs and Entrez Gene IDs

        Returns:
            pd.DataFrame: columns "drug_id" and "entrez_id"
        """
        df = pd.read_csv(f"{DataLoader.PATH_PDI}{DataLoader.PDI_DGIDB}", index_col=0)
        df["entrez_id"] = df["entrez_id"].map(DataLoader._clean_entrez)
        return df
