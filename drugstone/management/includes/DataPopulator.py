from drugstone.management.includes.DataLoader import DataLoader
import drugstone.models as models
from drugstone.management.includes.NodeCache import NodeCache


class DataPopulator:

    def __init__(self, cache: NodeCache):
        self.cache = cache

    def populate_expressions(self, update):
        if update:
            models.ExpressionLevel.objects.all().delete()

        self.cache.init_proteins()
        df = DataLoader.load_expressions()

        tissues_models = dict()
        for tissue_name in df.columns.values[2:]:
            try:
                tissue_model = models.Tissue.objects.get(name=tissue_name)
            except models.Tissue.DoesNotExist:
                tissue_model = models.Tissue.objects.create(name=tissue_name)
            tissues_models[tissue_name] = tissue_model

        proteins_linked = 0
        bulk = set()
        uniq = set()

        size = 0
        for _, row in df.iterrows():
            gene_name = row['Description']

            for protein_model in self.cache.get_proteins_by_gene(gene_name):
                proteins_linked += 1

                for tissue_name, tissue_model in tissues_models.items():
                    expr = models.ExpressionLevel(protein=protein_model,
                                                  tissue=tissue_model,
                                                  expression_level=row[tissue_name])
                    id = hash(expr)
                    if id in uniq:
                        continue
                    uniq.add(id)
                    bulk.add(expr)
            if len(bulk) > 100000:
                models.ExpressionLevel.objects.bulk_create(bulk)
                size += len(bulk)
                bulk = set()

        models.ExpressionLevel.objects.bulk_create(bulk)
        return size + len(bulk)

    def populate_ensg(self,update) -> int:
        """ Populates the Ensembl-Gene table in the django database.
        Also maps the added ensg entries to the corresponding proteins.
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many ensg-protein relations were added
        """
        if update:
            models.EnsemblGene.objects.all().delete()
        self.cache.init_proteins()
        data = DataLoader.load_ensg()
        bulk = list()

        for entrez, ensg_list in data.items():
            proteins = self.cache.get_proteins_by_entrez(entrez)
            for protein in proteins:
                for ensg in ensg_list:
                    bulk.append(models.EnsemblGene(name=ensg, protein=protein))
        models.EnsemblGene.objects.bulk_create(bulk)
        return len(bulk)

    def populate_ppi_string(self, dataset, update) -> int:
        """ Populates the Protein-Protein-Interactions from STRINGdb
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.cache.init_proteins()
        if update:
            models.ProteinProteinInteraction.objects.filter(ppi_dataset=dataset).delete()

        df = DataLoader.load_ppi_string()
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching proteins
                proteins_a = self.cache.get_proteins_by_entrez(row['entrez_a'])
                proteins_b = self.cache.get_proteins_by_entrez(row['entrez_b'])
            except KeyError:
                # continue if not found
                continue
            for protein_a in proteins_a:
                for protein_b in proteins_b:
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

    def populate_ppi_apid(self, dataset, update) -> int:
        """ Populates the Protein-Protein-Interactions from Apid
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.cache.init_proteins()

        if update:
            models.ProteinProteinInteraction.objects.filter(ppi_dataset=dataset).delete()
        df = DataLoader.load_ppi_apid()
        bulk = set()
        for _, row in df.iterrows():
            try:
                # try fetching proteins
                protein_a = self.cache.get_protein_by_uniprot(row['from_protein_ac'])
                protein_b = self.cache.get_protein_by_uniprot(row['to_protein_ac'])
            except KeyError:
                # continue if not found
                continue
            try:
                bulk.add(models.ProteinProteinInteraction(
                    ppi_dataset=dataset,
                    from_protein=protein_a,
                    to_protein=protein_b
                ))
            except models.ValidationError:
                continue
        models.ProteinProteinInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_ppi_biogrid(self,dataset, update) -> int:
        """ Populates the Protein-Protein-Interactions from BioGRID
        Handles loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.cache.init_proteins()

        if update:
            models.ProteinProteinInteraction.objects.filter(ppi_dataset=dataset).delete()
        df = DataLoader.load_ppi_biogrid()
        bulk = list()
        for _, row in df.iterrows():
            try:
                # try fetching proteins
                proteins_a = self.cache.get_proteins_by_entrez(row['entrez_a'])
                proteins_b = self.cache.get_proteins_by_entrez(row['entrez_b'])
            except KeyError:
                # TODO update error
                # continue if not found
                continue
            for protein_a in proteins_a:
                for protein_b in proteins_b:
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

    def populate_pdi_chembl(self,dataset, update) -> int:
        """ Populates the Protein-Drug-Interactions from Chembl
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.cache.init_proteins()
        self.cache.init_drugs()

        if update:
            models.ProteinDrugInteraction.objects.filter(pdi_dataset=dataset).delete()
        df = DataLoader.load_pdi_chembl()
        bulk = set()
        for _, row in df.iterrows():
            try:
                protein = self.cache.get_protein_by_uniprot(row['protein_ac'])
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = self.cache.get_drug_by_drugbank(row['drug_id'])
            except KeyError:
                # continue if not found
                continue
            bulk.add(models.ProteinDrugInteraction(
                pdi_dataset=dataset,
                protein=protein,
                drug=drug
            ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdis_disgenet(self, dataset, update) -> int:
        """ Populates the Protein-Disorder-Interactions from DisGeNET
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.cache.init_proteins()
        self.cache.init_disorders()


        if update:
            models.ProteinDisorderAssociation.objects.filter(pdis_dataset=dataset).delete()
        df = DataLoader.load_pdis_disgenet()
        bulk = set()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                protein = self.cache.get_protein_by_uniprot(row['protein_name'])
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching disorder
                disorder = self.cache.get_disorder_by_mondo(row['disorder_name'])
            except KeyError:
                # continue if not found
                continue
            bulk.add(models.ProteinDisorderAssociation(
                pdis_dataset=dataset,
                protein=protein,
                disorder=disorder,
                score=row['score']
            ))
        models.ProteinDisorderAssociation.objects.bulk_create(bulk)
        return len(bulk)

    def populate_drdis_drugbank(self, dataset, update) -> int:
        """ Populates the Drug-Disorder-Indications from DrugBank
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many edges were added
        """
        self.cache.init_drugs()
        self.cache.init_disorders()
        if update:
            models.DrugDisorderIndication.objects.filter(drdi_dataset=dataset).delete()

        df = DataLoader.load_drdis_drugbank()
        bulk = set()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                drug = self.cache.get_drug_by_drugbank(row['drugbank_id'])
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                disorder = self.cache.get_disorder_by_mondo(row['mondo_id'])
            except KeyError:
                # continue if not found
                continue
            bulk.add(models.DrugDisorderIndication(
                drdi_dataset=dataset,
                drug=drug,
                disorder=disorder,
            ))
        models.DrugDisorderIndication.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdi_dgidb(self,dataset, update) -> int:
        """ Populates the Protein-Drug-Interactions from DGIdb
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.cache.init_proteins()
        self.cache.init_drugs()

        if update:
            models.ProteinDrugInteraction.objects.filter(pdi_dataset=dataset).delete()
        df = DataLoader.load_pdi_dgidb()
        bulk = set()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                proteins = self.cache.get_proteins_by_entrez(row['entrez_id'])
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = self.cache.get_drug_by_drugbank(row['drug_id'])
            except KeyError:
                # continue if not found
                continue
            for protein in proteins:
                bulk.add(models.ProteinDrugInteraction(
                    pdi_dataset=dataset,
                    protein=protein,
                    drug=drug
                ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)

    def populate_pdi_drugbank(self,dataset, update) -> int:
        """ Populates the Protein-Drug-Interactions from Drugbank
        Handles Loading the data and passing it to the django database

        Returns:
            int: Count of how many interactions were added
        """
        self.cache.init_proteins()
        self.cache.init_drugs()


        if update:
            models.ProteinDrugInteraction.objects.filter(pdi_dataset=dataset).delete()
        df = DataLoader.load_pdi_drugbank()
        bulk = set()
        for _, row in df.iterrows():
            try:
                # try fetching protein
                proteins = self.cache.get_proteins_by_entrez(row['entrez_id'])
            except KeyError:
                # continue if not found
                continue
            try:
                # try fetching drug
                drug = self.cache.get_drug_by_drugbank(row['drug_id'])
            except KeyError:
                # continue if not found
                continue
            for protein in proteins:
                bulk.add(models.ProteinDrugInteraction(
                    pdi_dataset=dataset,
                    protein=protein,
                    drug=drug
                ))
        models.ProteinDrugInteraction.objects.bulk_create(bulk)
        return len(bulk)
