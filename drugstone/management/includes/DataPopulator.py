from drugstone.management.includes.DataLoader import DataLoader
import drugstone.models as models


class DataPopulator:

    def populate_proteins() -> int:
        """ Populates the Protein table in the django database.
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many proteins were added
        """
        df = DataLoader.load_proteins()
        proteins = list()
        for _, row in df.iterrows():
            proteins.append(models.Protein(
                uniprot_code=row['protein_ac'],
                gene=row['gene_name'],
                entrez=row['entrez_id'],
                protein_name=row['protein_name'])
            )
        models.Protein.objects.bulk_create(proteins)
        return len(proteins)

    def populate_disorders() -> int:
        """ Populates the Disorder table in the django database.
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many disorders were added
        """
        df = DataLoader.load_disorders()
        bulk = list()
        for _, row in df.iterrows():
            bulk.append(models.Disorder(
                mondo_id=row['mondo_id'],
                label=row['label'],
                icd10=row['icd10']
            ))
        models.Disorder.objects.bulk_create(bulk)
        return len(bulk)

    def populate_ensg() -> int:
        """ Populates the Ensembl-Gene table in the django database.
        Also maps the added ensg entries to the corresponding proteins.
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many ensg-protein relations were added
        """
        data = DataLoader.load_ensg()
        bulk = list()
        for entrez, ensg_list in data.items():
            protein = models.Protein.objects.get(entrez=entrez)
            for ensg in ensg_list:
                bulk.append(models.EnsemblGene(name=ensg, protein=protein))
        models.EnsemblGene.objects.bulk_create(bulk)
        return len(bulk)

    def populate_ppi_string() -> int:
        """ Populates the Protein-Protein-Interactions from STRINGdb
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
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
                protein_a = models.Protein.objects.get(entrez=row['entrez_a'])
                protein_b = models.Protein.objects.get(entrez=row['entrez_b'])
            except models.Protein.DoesNotExist:
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

    def populate_ppi_apid() -> int:
        """ Populates the Protein-Protein-Interactions from Apid
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
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
                protein_a = models.Protein.objects.get(uniprot_code=row['from_protein_ac'])
                protein_b = models.Protein.objects.get(uniprot_code=row['to_protein_ac'])
            except models.Protein.DoesNotExist:
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

    def populate_ppi_biogrid() -> int:
        """ Populates the Protein-Protein-Interactions from BioGRID
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
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
                protein_a = models.Protein.objects.get(entrez=row['entrez_a'])
                protein_b = models.Protein.objects.get(entrez=row['entrez_b'])
            except models.Protein.DoesNotExist:
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

    def populate_pdi_chembl() -> int:
        """ Populates the Protein-Drug-Interactions from Chembl
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        df = DataLoader.load_pdi_chembl()
        dataset, _ = models.PDIDataset.objects.get_or_create(
            name='ChEMBL',
            link='https://www.ebi.ac.uk/chembl/',
            version='27',
        )
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                protein = models.Protein.objects.get(uniprot_code=row['protein_ac'])
            except models.Protein.DoesNotExist:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = models.Drug.objects.get(drug_id=row['drug_id'])
            except models.Drug.DoesNotExist:
                # continue if not found
                continue
            bulk.append(models.ProteinDrugInteraction(
                pdi_dataset=dataset,
                protein=protein,
                drug=drug
            ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdis_disgenet() -> int:
        """ Populates the Protein-Disorder-Interactions from DisGeNET
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
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
                protein = models.Protein.objects.get(uniprot_code=row['protein_name'])
            except models.Protein.DoesNotExist:
                # continue if not found
                continue
            try:
                # try fetching drug
                disorder = models.Disorder.objects.get(mondo_id=row['disorder_name'])
            except models.Disorder.DoesNotExist:
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

    def populate_drdis_drugbank() -> int:
        """ Populates the Drug-Disorder-Indications from DrugBank
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many edges were added
        """
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
                drug = models.Drug.objects.get(drug_id=row['drugbank_id'])
            except models.Drug.DoesNotExist:
                # continue if not found
                continue
            try:
                # try fetching drug
                disorder = models.Disorder.objects.get(mondo_id=row['mondo_id'])
            except models.Disorder.DoesNotExist:
                # continue if not found
                continue
            bulk.append(models.DrugDisorderIndication(
                drdi_dataset=dataset,
                drug=drug,
                disorder=disorder,
            ))
        models.DrugDisorderIndication.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdi_dgidb() -> int:
        """ Populates the Protein-Drug-Interactions from DGIdb
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
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
                protein = models.Protein.objects.get(entrez=row['entrez_id'])
            except models.Protein.DoesNotExist:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = models.Drug.objects.get(drug_id=row['drug_id'])
            except models.Drug.DoesNotExist:
                # continue if not found
                continue
            bulk.append(models.ProteinDrugInteraction(
                pdi_dataset=dataset,
                protein=protein,
                drug=drug
            ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdi_drugbank() -> int:
        """ Populates the Protein-Drug-Interactions from Drugbank
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
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
                protein = models.Protein.objects.get(entrez=row['entrez_id'])
            except models.Protein.DoesNotExist:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = models.Drug.objects.get(drug_id=row['drug_id'])
            except models.Drug.DoesNotExist:
                # continue if not found
                continue
            bulk.append(models.ProteinDrugInteraction(
                pdi_dataset=dataset,
                protein=protein,
                drug=drug
            ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)
