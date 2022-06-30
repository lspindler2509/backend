from collections import defaultdict
import drugstone.models as models


class NodeCache:

    proteins = dict()
    entrez_to_uniprot = defaultdict(lambda: set())
    gene_name_to_uniprot = defaultdict(lambda: set())
    disorders = dict()
    drugs = dict()

    def init_protein_maps(self):
        print("Generating protein id maps...")
        for protein in self.proteins.values():
            self.entrez_to_uniprot[protein.entrez].add(protein.uniprot_code)
            self.gene_name_to_uniprot[protein.gene].add(protein.uniprot_code)

    def init_proteins(self):
        if len(self.proteins) == 0:
            print("Generating protein maps...")
            for protein in models.Protein.objects.all():
                self.proteins[protein.uniprot_code] = protein
        if len(self.proteins) > 0 and (len(self.entrez_to_uniprot) == 0 or len(self.gene_name_to_uniprot) == 0):
            self.init_protein_maps()


    def init_drugs(self):
        if len(self.drugs) == 0:
            print("Generating drug map...")
            for drug in models.Drug.objects.all():
                self.drugs[drug.drug_id] = drug

    def init_disorders(self):
        if len(self.disorders) == 0:
            print("Generating disorder map...")
            for disorder in models.Disorder.objects.all():
                self.disorders[disorder.mondo_id] = disorder


    def get_protein_by_uniprot(self,uniprot_id):
        return self.proteins[uniprot_id]

    def get_proteins_by_entrez(self,entrez_id):
        out = list()
        for g in self.entrez_to_uniprot[entrez_id]:
            out.append(self.proteins[g])
        return out

    def get_proteins_by_gene(self, gene_name):
        out = list()
        for g in self.gene_name_to_uniprot[gene_name]:
            out.append(self.proteins[g])
        return out

    def get_drug_by_drugbank(self, drugbank_id):
        return self.drugs[drugbank_id]

    def get_disorder_by_mondo(self, mondo_id):
        return self.disorders[mondo_id]