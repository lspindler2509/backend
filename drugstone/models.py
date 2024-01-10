from django.core.exceptions import ValidationError
from django.db import models


# Main biological and medical entities


class PPIDataset(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=False)
    link = models.CharField(max_length=128, default="", unique=False)
    version = models.CharField(max_length=128, default="", unique=False)
    licenced = models.BooleanField(default=True)

    def __str__(self):
        return f'{self.name}-{self.version}_{"licenced" if self.licenced else "unlicenced"}'
    
    class Meta:
        unique_together = ("name", "version", "licenced")


class PDIDataset(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=False)
    link = models.CharField(max_length=128, default="", unique=False)
    version = models.CharField(max_length=128, default="", unique=False)
    licenced = models.BooleanField(default=True)

    def __str__(self):
        return f'{self.name}-{self.version}_{"licenced" if self.licenced else "unlicenced"}'

    class Meta:
        unique_together = ("name", "version", "licenced")


class PDisDataset(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=False)
    link = models.CharField(max_length=128, default="", unique=False)
    version = models.CharField(max_length=128, default="", unique=False)
    licenced = models.BooleanField(default=True)

    def __str__(self):
        return f'{self.name}-{self.version}_{"licenced" if self.licenced else "unlicenced"}'

    class Meta:
        unique_together = ("name", "version", "licenced")


class DrDiDataset(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=False)
    link = models.CharField(max_length=128, default="", unique=False)
    version = models.CharField(max_length=128, default="", unique=False)
    licenced = models.BooleanField(default=True)

    def __str__(self):
        return f'{self.name}-{self.version}_{"licenced" if self.licenced else "unlicenced"}'

    class Meta:
        unique_together = ("name", "version", "licenced")


class EnsemblGene(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=15)  # starts with ENSG...
    protein = models.ForeignKey(
        "Protein", on_delete=models.CASCADE, related_name="ensg"
    )
    

class Protein(models.Model):
    # According to https://www.uniprot.org/help/accession_numbers UniProt accession codes
    # are either 6 or 10 characters long
    id = models.AutoField(primary_key=True)
    uniprot_code = models.CharField(max_length=10)
    gene = models.CharField(max_length=127, default="")  # symbol
    protein_name = models.CharField(max_length=255, default="")
    entrez = models.CharField(max_length=15, default="")
    drugs = models.ManyToManyField(
        "Drug", through="ProteinDrugInteraction", related_name="interacting_drugs"
    )
    tissue_expression = models.ManyToManyField(
        "Tissue", through="ExpressionLevel", related_name="interacting_drugs"
    )

    class Meta:
        unique_together = ("uniprot_code", "gene", "entrez")

    def __str__(self):
        return self.gene

    def __eq__(self, other):
        return (
            self.uniprot_code == other.uniprot_code
            and self.gene == other.gene
            and self.protein_name == other.protein_name
            and self.entrez == other.entrez
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.uniprot_code, self.gene, self.entrez))

    def update(self, other):
        self.uniprot_code = other.uniprot_code
        self.gene = other.gene
        self.protein_name = other.protein_name
        self.entrez = other.entrez
        

class ExpressionLevel(models.Model):
    id = models.AutoField(primary_key=True)
    tissue = models.ForeignKey("Tissue", on_delete=models.CASCADE)
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE)
    expression_level = models.FloatField()

    class Meta:
        unique_together = ("tissue", "protein")

    def __hash__(self):
        return hash(f"{self.tissue_id}_{self.protein_id}")


class Tissue(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=True)

    def __str__(self):
        return self.name


class Disorder(models.Model):
    id = models.AutoField(primary_key=True)
    mondo_id = models.CharField(max_length=7)
    label = models.CharField(max_length=256, default="")  # symbol
    icd10 = models.CharField(max_length=512, default="")
    proteins = models.ManyToManyField(
        "Protein",
        through="ProteinDisorderAssociation",
        related_name="associated_proteins",
    )

    class Meta:
        unique_together = ("mondo_id", "label", "icd10")

    def __str__(self):
        return self.label

    def __eq__(self, other):
        return (
            self.mondo_id == other.mondo_id
            and self.label == other.label
            and self.icd10 == other.icd10
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.mondo_id, self.label, self.icd10))

    def update(self, other):
        self.mondo_id = other.mondo_id
        self.label = other.label
        self.icd10 = other.icd10
        
        
class Drug(models.Model):
    id = models.AutoField(primary_key=True)
    drug_id = models.CharField(max_length=10, unique=True)
    name = models.CharField(max_length=256, default="")
    status = models.CharField(max_length=128, default="")
    # in_trial = models.BooleanField(default=False)
    # in_literature = models.BooleanField(default=False)
    links = models.CharField(max_length=16 * 1024, default="")

    def __str__(self):
        return self.drug_id

    def __eq__(self, other):
        return (
            self.drug_id == other.drug_id
            and self.name == other.name
            and self.status == other.status
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.drug_id, self.name, self.status))

    def update(self, other):
        self.drug_id = other.drug_id
        self.name = other.name
        self.status = other.status
        self.links = other.links


class ProteinDisorderAssociationAbstract(models.Model):
    drugstone_id = models.BigAutoField(primary_key=True)
    pdis_dataset = models.ForeignKey(
        "PDisDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="protein_disorder_associations"
    )
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE, related_name="disorder_associations")
    disorder = models.ForeignKey("Disorder", on_delete=models.CASCADE, related_name="protein_associations")
    score = models.FloatField()

    class Meta:
        unique_together = ("pdis_dataset", "protein", "disorder")
        ordering = ('pdis_dataset__name', '-pdis_dataset__id')
        abstract = True

    def __str__(self):
        return f"{self.pdis_dataset}-{self.protein}-{self.disorder}"

    def __eq__(self, other):
        return (
            self.pdis_dataset.name == other.pdis_dataset.name
            and self.pdis_dataset.licenced == other.pdis_dataset.licenced
            and self.protein == other.protein
            and self.disorder == other.disorder
            and self.score == self.score
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.pdis_dataset.name, self.pdis_dataset.licenced, self.protein, self.disorder, self.score))


class ProteinDisorderAssociationManager(models.Manager):
    def bulk_create(self, objs: set, removed: set = set(), **kwargs):
        if len(objs) or len(removed):
            
            history_objs = [*create_history_objs_new(objs, ProteinDisorderAssociationHistory), 
                    *create_history_objs_removed(removed, ProteinDisorderAssociationHistory, 'pdis_dataset', PDIDataset)]
            # store changes in history
            ProteinDisorderAssociationHistory.objects.bulk_create(history_objs, **kwargs)
        return super().bulk_create(objs, **kwargs)
        
        
class ProteinDisorderAssociation(ProteinDisorderAssociationAbstract):
    objects = ProteinDisorderAssociationManager()


class DrugDisorderIndicationAbstract(models.Model):
    drugstone_id = models.AutoField(primary_key=True)
    drdi_dataset = models.ForeignKey(
        "DrDiDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="drug_disorder_indications"
    )
    drug = models.ForeignKey("Drug", on_delete=models.CASCADE, related_name="disorder_indications")
    disorder = models.ForeignKey("Disorder", on_delete=models.CASCADE, related_name="drug_indications")

    class Meta:
        unique_together = ("drdi_dataset", "drug", "disorder")
        abstract = True
        ordering = ('drdi_dataset__name', '-drdi_dataset__id')
        
    def __str__(self):
        return f"{self.drdi_dataset}-{self.drug}-{self.disorder}"

    def __eq__(self, other):
        return (
            self.drdi_dataset.name == other.drdi_dataset.name
            and self.drdi_dataset.licenced == other.drdi_dataset.licenced
            and self.drug == other.drug
            and self.disorder == other.disorder
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.drdi_dataset.name, self.drdi_dataset.licenced, self.drug, self.disorder))


class DrugDisorderIndicationManager(models.Manager):
    def bulk_create(self, objs: set, removed: set = set(), **kwargs):
        if len(objs) or len(removed):
            history_objs = [*create_history_objs_new(objs, DrugDisorderIndicationHistory), 
                *create_history_objs_removed(removed, DrugDisorderIndicationHistory, 'drdi_dataset', DrDiDataset)]
            # store changes in history
            DrugDisorderIndicationHistory.objects.bulk_create(history_objs, **kwargs)
        return super().bulk_create(objs, **kwargs)


class DrugDisorderIndication(DrugDisorderIndicationAbstract):
    objects = DrugDisorderIndicationManager()


class ProteinProteinInteractionAbstract(models.Model):
    drugstone_id = models.BigAutoField(primary_key=True)
    ppi_dataset = models.ForeignKey(
        "PPIDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="protein_interations"
    )
    from_protein = models.ForeignKey(
        "Protein", on_delete=models.CASCADE, 
        related_name="protein_interactions_out"
    )
    to_protein = models.ForeignKey(
        "Protein", on_delete=models.CASCADE, 
        related_name="protein_interactions_to"
    )

    def __str__(self):
        return f"{self.ppi_dataset}-{self.from_protein}-{self.to_protein}"

    def __eq__(self, other):
        return (
            self.ppi_dataset.name == other.ppi_dataset.name
            and self.ppi_dataset.licenced == other.ppi_dataset.licenced
            and self.from_protein == other.from_protein
            and self.to_protein == other.to_protein
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.ppi_dataset.name, self.ppi_dataset.licenced, self.from_protein, self.to_protein))
    
    class Meta:
        abstract = True
        ordering = ('ppi_dataset__name', '-ppi_dataset__id')


class ProteinProteinInteractionManager(models.Manager):
    def bulk_create(self, objs: set, removed: set = set(), **kwargs):
        if len(objs) or len(removed):
            history_objs = [*create_history_objs_new(objs, ProteinProteinInteractionHistory), 
                *create_history_objs_removed(removed, ProteinProteinInteractionHistory, 'ppi_dataset', PPIDataset)]
            # store changes in history
            ProteinProteinInteractionHistory.objects.bulk_create(history_objs, **kwargs)
        return super().bulk_create(objs, **kwargs)    

    
class ProteinProteinInteraction(ProteinProteinInteractionAbstract):
    objects = ProteinProteinInteractionManager()
    
    def validate_unique(self, exclude=None):
        p1p2_q = ProteinProteinInteraction.objects.filter(
            from_protein=self.from_protein, 
            to_protein=self.to_protein,
            ppi_dataset=self.ppi_dataset,
        )
        p2p1_q = ProteinProteinInteraction.objects.filter(
            from_protein=self.to_protein,
            to_protein=self.from_protein,
            ppi_dataset=self.ppi_dataset,
        )

        if p1p2_q.exists() or p2p1_q.exists():
            raise ValidationError("Protein-Protein interaction must be unique!")
    
    def save(self, *args, **kwargs):
        self.validate_unique()
        super(ProteinProteinInteraction, self).save(*args, **kwargs)


class ProteinDrugInteractionAbstract(models.Model):
    drugstone_id = models.BigAutoField(primary_key=True)
    pdi_dataset = models.ForeignKey(
        PDIDataset,
        null=True,
        on_delete=models.CASCADE,
        related_name="protein_drug_interations"
    )
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE, related_name="drug_interaction")
    drug = models.ForeignKey("Drug", on_delete=models.CASCADE, related_name="drug_interaction")
    actions = models.CharField(max_length=255, default="[]")

    class Meta:
        unique_together = ("pdi_dataset", "protein", "drug")
        ordering = ('pdi_dataset__name', '-pdi_dataset__id')
        abstract = True

    def __str__(self):
        return f"{self.pdi_dataset}-{self.protein}-{self.drug}"

    def __eq__(self, other):
        return (
            self.pdi_dataset.name == other.pdi_dataset.name
            and self.pdi_dataset.licenced == other.pdi_dataset.licenced
            and self.protein == other.protein
            and self.drug == other.drug
            and self.actions == other.actions
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.pdi_dataset.name, self.pdi_dataset.licenced, self.protein, self.drug, self.actions))


class ProteinDrugInteractionManager(models.Manager):
    def bulk_create(self, objs: set, removed: set = set(), **kwargs):
        bulk_create_response = super().bulk_create(objs, **kwargs)
        if len(objs) or len(removed):
            history_objs = [*create_history_objs_new(objs, ProteinDrugInteractionHistory), 
                *create_history_objs_removed(removed, ProteinDrugInteractionHistory, 'pdi_dataset', PDIDataset)]
            # store changes in history
            ProteinDrugInteractionHistory.objects.bulk_create(history_objs, **kwargs)
        return bulk_create_response
    

class ProteinDrugInteraction(ProteinDrugInteractionAbstract):
    objects = ProteinDrugInteractionManager()
    
    def historic(self, objs, historic_dataset: PPIDataset):
        # load historic version of objs
        historic_changes = ProteinDrugInteractionHistory.filter(pdi_dataset__name=historic_dataset.name)
        # convert list of objs to dict for faster lookup, key is composed out of unique identifier for element
        objs = left_search(objs, historic_changes, historic_dataset)
        return objs


def left_search(objs, historic_changes, historic_dataset, key_attributes):
    objs_dict = {hash((getattr(change, attr) for attr in key_attributes)): obj for obj in objs}
    left_search_done = set()
    # changes are sorted by version, hence in correct order
    for change in historic_changes:
        key = hash((getattr(change, attr) for attr in key_attributes))
        if key in left_search_done:
            # left-search is done for the element affected by this change
            continue
        
        if change.new and change.pdi_dataset.id > historic_dataset.id:
            # item was added in this version, hence it did not exist in previous version
            del objs_dict[key]
        elif change.pdi_dataset.id <= historic_dataset.id:
            # update and return, left-search terminated, state of target dataset reached
            objs_dict[key] = change
            left_search_done.add(key)
        else:
            # change.new --> False && change.pdi_dataset.id <= historic_dataset.id
            # update but continue left search
            objs_dict[key] = change
    return objs_dict.values()


class Task(models.Model):
    token = models.CharField(max_length=32, unique=True, primary_key=True)
    created_at = models.DateTimeField(auto_now_add=True)
    target = models.CharField(
        max_length=32, choices=[("drug", "Drug"), ("drug-target", "Drug Target")]
    )

    algorithm = models.CharField(max_length=128)
    parameters = models.TextField()

    progress = models.FloatField(default=0.0)  # Progress as fraction (0.0 - 1.0)
    started_at = models.DateTimeField(null=True)
    finished_at = models.DateTimeField(null=True)
    worker_id = models.CharField(max_length=128, null=True)
    job_id = models.CharField(max_length=128, null=True)
    done = models.BooleanField(default=False)
    failed = models.BooleanField(default=False)
    status = models.CharField(max_length=255, null=True)

    result = models.TextField(null=True)


class Network(models.Model):
    id = models.CharField(primary_key=True, max_length=32, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    nodes = models.TextField(null=True, default="")
    edges = models.TextField(null=True, default="")
    config = models.TextField(null=True, default="")
    groups = models.TextField(null=True, default="")


def create_history_objs_new(objs, history_model):
    history_objs_changed = []
    for hobj in objs:
        history_model_fields = [f.name for f in history_model._meta.fields]
        history_model_dict = {x: getattr(hobj, x, None) for x in history_model_fields}
        hobj = history_model(**history_model_dict)
        hobj.new = True
        hobj.removed = False
        history_objs_changed.append(hobj)
    return history_objs_changed
    
    
def create_history_objs_removed(objs_removed, history_model, dataset_key: str, dataset_model):
    history_objs_changed = []
    dataset_map = {}
    for hobj in objs_removed:
        history_model_fields = [f.name for f in history_model._meta.fields]
        history_model_dict = {x: getattr(hobj, x, None) for x in history_model_fields}
        hobj = history_model(**history_model_dict)
        hobj.new = False
        hobj.removed = True
        
        # update to latets dataset version for history table
        dataset_old = getattr(hobj, dataset_key)
        if str(dataset_old) in dataset_map:
            dataset_new = dataset_map[str(dataset_old)]
        else:
            dataset_new = dataset_model.objects.filter(name=dataset_old.name, licenced=dataset_old.licenced).order_by('-id')[0]
            # save to avoid unnecessary db lookups
            dataset_map[str(dataset_old)] = dataset_new
        setattr(hobj, dataset_key, dataset_new)
        
        history_objs_changed.append(hobj)
    return history_objs_changed

    
class ProteinProteinInteractionHistory(ProteinProteinInteractionAbstract):
    id = models.AutoField(primary_key=True)
    drugstone_id = models.BigIntegerField()
    new = models.BooleanField(default=True)
    removed = models.BooleanField(default=False)
    
    # to set different related name
    ppi_dataset = models.ForeignKey(
        "PPIDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="historic_ppi_interactions"
    )
    from_protein = models.ForeignKey(
        "Protein", on_delete=models.CASCADE, 
        related_name="historic_protein_interactions_out"
    )
    to_protein = models.ForeignKey(
        "Protein", on_delete=models.CASCADE, 
        related_name="historic_protein_interactions_in"
    )

    
class ProteinDisorderAssociationHistory(ProteinDisorderAssociationAbstract):
    id = models.AutoField(primary_key=True)
    pdis_dataset = models.ForeignKey(
        "PDisDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="historic_protein_disorder_associations"
    )
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE, related_name="historic_disorder_associations")
    disorder = models.ForeignKey("Disorder", on_delete=models.CASCADE, related_name="historic_protein_associations")
    drugstone_id = models.BigIntegerField()
    new = models.BooleanField(default=True)
    removed = models.BooleanField(default=False)   
    
    
class ProteinDrugInteractionHistory(ProteinDrugInteractionAbstract):
    id = models.AutoField(primary_key=True)
    pdi_dataset = models.ForeignKey(
        PDIDataset,
        null=True,
        on_delete=models.CASCADE,
        related_name="historic_protein_drug_interations"
    )
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE, related_name="historic_drug_interaction")
    drug = models.ForeignKey("Drug", on_delete=models.CASCADE, related_name="historic_drug_interaction")
    drugstone_id = models.BigIntegerField()
    new = models.BooleanField(default=True)
    removed = models.BooleanField(default=False)


class DrugDisorderIndicationHistory(DrugDisorderIndicationAbstract):
    id = models.AutoField(primary_key=True)
    drdi_dataset = models.ForeignKey(
        "DrDiDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="historic_drug_disorder_indications"
    )
    drug = models.ForeignKey("Drug", on_delete=models.CASCADE, related_name="historic_disorder_indications")
    disorder = models.ForeignKey("Disorder", on_delete=models.CASCADE, related_name="historic_drug_indications")
    drugstone_id = models.BigIntegerField()
    new = models.BooleanField(default=True)
    removed = models.BooleanField(default=False)
