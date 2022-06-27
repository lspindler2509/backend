from collections import defaultdict

from drugstone.management.includes.DataLoader import DataLoader
import drugstone.models as models


class DataPopulator:
    proteins = dict()
    uniprot_to_ensembl = dict()
    gene_name_to_ensembl = defaultdict(lambda: set())
    disorders = dict()
    drugs = dict()

    def init_proteins(self):
        if len(self.proteins) == 0:
            print("Generating protein maps...")
            for protein in models.Protein.objects.all():
                self.proteins[protein.entrez]=protein
                self.uniprot_to_ensembl[protein.uniprot_code] = protein.entrez
                self.gene_name_to_ensembl[protein.gene].add(protein.entrez)

    def init_drugs(self):
        if len(self.drugs)== 0:
            print("Generating drug map...")
            for drug in models.Drug.objects.all():
                self.drugs[drug.drug_id]=drug

    def init_disorders(self):
        if len(self.disorders) == 0:
            print("Generating disorder map...")
            for disorder in models.Disorder.objects.all():
                self.disorders[disorder.mondo_id]=disorder

    # def populate_proteins(self) -> int:
    #     """ Populates the Protein table in the django database.
    #     Handles loading the data and passing it to the django database
    #
    #     Returns:
    #         int: Count of how many proteins were added
    #     """
    #     df = DataLoader.load_proteins()
    #     for _, row in df.iterrows():
    #         self.proteins[row['entrez_id']] = models.Protein(
    #             uniprot_code=row['protein_ac'],
    #             gene=row['gene_name'],
    #             entrez=row['entrez_id'],
    #             protein_name=row['protein_name'])
    #         self.uniprot_to_ensembl[row['protein_ac']] = row['entrez_id']
    #         self.gene_name_to_ensembl[row['gene_name']].add(row['entrez_id'])
    #
    #     models.Protein.objects.bulk_create(self.proteins.values())
    #     return len(self.proteins)
    #
    # def populate_disorders(self) -> int:
    #     """ Populates the Disorder table in the django database.
    #     Handles loading the data and passing it to the django database
    #
    #     Returns:
    #         int: Count of how many disorders were added
    #     """
    #     df = DataLoader.load_disorders()
    #     for _, row in df.iterrows():
    #         self.disorders[row['mondo_id']] = models.Disorder(
    #             mondo_id=row['mondo_id'],
    #             label=row['label'],
    #             icd10=row['icd10']
    #         )
    #     models.Disorder.objects.bulk_create(self.disorders.values())
    #     return len(self.disorders)
    #
    # def populate_drugs(self):
    #     df = DataLoader.load_drugs()
    #     for _, row in df.iterrows():
    #         drug_id = row['drug_id']
    #         drug_name = row['drug_name']
    #         drug_status = row['drug_status']
    #         self.drugs[drug_id] = models.Drug(
    #             drug_id=drug_id,
    #             name=drug_name,
    #             status=drug_status)
    #     models.Drug.objects.bulk_create(self.drugs.values())
    #     return len(self.drugs)

    def populate_expessions(self):
        self.init_proteins()
        df = DataLoader.load_expressions()

        tissues_models = dict()
        for tissue_name in df.columns.values[2:]:
            try:
                tissue_model = models.Tissue.objects.get(name=tissue_name)
            except models.Tissue.DoesNotExist:
                tissue_model = models.Tissue.objects.create(name=tissue_name)
            tissues_models[tissue_name] = tissue_model

        proteins_linked = 0
        unique = set()
        bulk = list()

        for _, row in df.iterrows():
            gene_name = row['Description']

            for protein_id in self.gene_name_to_ensembl[gene_name]:
                protein_model = self.proteins[protein_id]
                proteins_linked += 1

                for tissue_name, tissue_model in tissues_models.items():
                    id = f"{tissue_name}_{protein_id}"
                    if id in unique:
                        continue
                    unique.add(id)
                    bulk.append(models.ExpressionLevel(protein=protein_model,
                                                       tissue=tissue_model,
                                                       expression_level=row[tissue_name]))
        models.ExpressionLevel.objects.bulk_create(bulk)
        return len(bulk)

    def populate_ensg(self) -> int:
        """ Populates the Ensembl-Gene table in the django database.
        Also maps the added ensg entries to the corresponding proteins.
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many ensg-protein relations were added
        """
        self.init_proteins()
        data = DataLoader.load_ensg()
        bulk = list()
        for entrez, ensg_list in data.items():
            protein = self.proteins[entrez]
            for ensg in ensg_list:
                bulk.append(models.EnsemblGene(name=ensg, protein=protein))
        models.EnsemblGene.objects.bulk_create(bulk)
        return len(bulk)

    def populate_ppi_string(self) -> int:
        """ Populates the Protein-Protein-Interactions from STRINGdb
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.init_proteins()
        df = DataLoader.load_ppi_string()
        dataset, _ = models.PPIDataset.objects.get_or_create(
            name='STRING',
            link='https://string-db.org/',
            version='11.0'
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching proteins
                protein_a = self.proteins[row['entrez_a']]
                protein_b = self.proteins[row['entrez_b']]
            except KeyError:
                # continue if not found
                continue
            try:
                bulk.append(models.ProteinProteinInteraction(
                    ppi_dataset=dataset,
                    from_protein=protein_a,
                    to_protein=protein_b
                ))
            except models.ValidationError:
                # duplicate
                continue
        models.ProteinProteinInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_ppi_apid(self) -> int:
        """ Populates the Protein-Protein-Interactions from Apid
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.init_proteins()
        df = DataLoader.load_ppi_apid()
        dataset, _ = models.PPIDataset.objects.get_or_create(
            name='APID',
            link='http://cicblade.dep.usal.es:8080/APID/',
            version='January 2019'
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching proteins
                protein_a = self.proteins[self.uniprot_to_ensembl[row['from_protein_ac']]]
                protein_b = self.proteins[self.uniprot_to_ensembl[row['to_protein_ac']]]
            except KeyError:
                # continue if not found
                continue
            try:
                bulk.append(models.ProteinProteinInteraction(
                    ppi_dataset=dataset,
                    from_protein=protein_a,
                    to_protein=protein_b
                ))
            except models.ValidationError:
                continue
        models.ProteinProteinInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_ppi_biogrid(self) -> int:
        """ Populates the Protein-Protein-Interactions from BioGRID
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.init_proteins()
        df = DataLoader.load_ppi_biogrid()
        dataset, _ = models.PPIDataset.objects.get_or_create(
            name='BioGRID',
            link='https://thebiogrid.org/',
            version='4.0'
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching proteins
                protein_a = self.proteins[row['entrez_a']]
                protein_b = self.proteins[row['entrez_b']]
            except KeyError:
                # TODO update error
                # continue if not found
                continue
            try:
                bulk.append(models.ProteinProteinInteraction(
                    ppi_dataset=dataset,
                    from_protein=protein_a,
                    to_protein=protein_b
                ))
            except models.ValidationError:
                # duplicate
                continue
        models.ProteinProteinInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdi_chembl(self) -> int:
        """ Populates the Protein-Drug-Interactions from Chembl
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.init_proteins()
        self.init_drugs()
        df = DataLoader.load_pdi_chembl()
        dataset, _ = models.PDIDataset.objects.get_or_create(
            name='ChEMBL',
            link='https://www.ebi.ac.uk/chembl/',
            version='27',
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                protein = self.proteins[self.uniprot_to_ensembl[row['protein_ac']]]
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = self.drugs[row['drug_id']]
            except KeyError:
                # continue if not found
                continue
            bulk.append(models.ProteinDrugInteraction(
                pdi_dataset=dataset,
                protein=protein,
                drug=drug
            ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdis_disgenet(self,) -> int:
        """ Populates the Protein-Disorder-Interactions from DisGeNET
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.init_proteins()
        self.init_disorders()
        df = DataLoader.load_pdis_disgenet()
        dataset, _ = models.PDisDataset.objects.get_or_create(
            name='DisGeNET',
            link='https://www.disgenet.org/home/',
            version='6.0',
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                protein = self.proteins[self.uniprot_to_ensembl[row['protein_name']]]
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                disorder = self.disorders[str(int(row['disorder_name']))]
            except KeyError:
                # continue if not found
                continue
            bulk.append(models.ProteinDisorderAssociation(
                pdis_dataset=dataset,
                protein=protein,
                disorder=disorder,
                score=row['score']
            ))
        models.ProteinDisorderAssociation.objects.bulk_create(bulk)
        return len(bulk)

    def populate_drdis_drugbank(self) -> int:
        """ Populates the Drug-Disorder-Indications from DrugBank
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many edges were added
        """
        self.init_drugs()
        self.init_disorders()
        df = DataLoader.load_drdis_drugbank()
        dataset, _ = models.DrDiDataset.objects.get_or_create(
            name='DrugBank',
            link='https://go.drugbank.com/',
            version='5.1.8',
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                drug = self.drugs[row['drugbank_id']]
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                disorder = self.disorders[str(int(row['mondo_id']))]
            except KeyError:
                # continue if not found
                continue
            bulk.append(models.DrugDisorderIndication(
                drdi_dataset=dataset,
                drug=drug,
                disorder=disorder,
            ))
        models.DrugDisorderIndication.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdi_dgidb(self) -> int:
        """ Populates the Protein-Drug-Interactions from DGIdb
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.init_proteins()
        self.init_drugs()
        df = DataLoader.load_pdi_dgidb()
        dataset, _ = models.PDIDataset.objects.get_or_create(
            name='DGIdb',
            link='https://www.dgidb.org/',
            version='4.2.0'
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                protein = self.proteins[row['entrez_id']]
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = self.drugs[row['drug_id']]
            except KeyError:
                # continue if not found
                continue
            bulk.append(models.ProteinDrugInteraction(
                pdi_dataset=dataset,
                protein=protein,
                drug=drug
            ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdi_drugbank(self) -> int:
        """ Populates the Protein-Drug-Interactions from Drugbank
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.init_proteins()
        self.init_drugs()
        df = DataLoader.load_pdi_drugbank()
        dataset, _ = models.PDIDataset.objects.get_or_create(
            name='DrugBank',
            link='https://go.drugbank.com/',
            version='5.1.7'
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                protein = self.proteins[row['entrez_id']]
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = self.drugs[row['drug_id']]
            except KeyError:
                # continue if not found
                continue
            bulk.append(models.ProteinDrugInteraction(
                pdi_dataset=dataset,
                protein=protein,
                drug=drug
            ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)
