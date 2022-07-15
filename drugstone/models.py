from django.core.exceptions import ValidationError
from django.db import models


# Main biological and medical entities


class Tissue(models.Model):
    name = models.CharField(max_length=128, default='', unique=True)

    def __str__(self):
        return self.name


class PPIDataset(models.Model):
    name = models.CharField(max_length=128, default='', unique=False)
    link = models.CharField(max_length=128, default='', unique=False)
    version = models.CharField(max_length=128, default='', unique=False)

    def __str__(self):
        return f'{self.name}-{self.version}'

    class Meta:
        unique_together = ('name', 'version')


class PDIDataset(models.Model):
    name = models.CharField(max_length=128, default='', unique=False)
    link = models.CharField(max_length=128, default='', unique=False)
    version = models.CharField(max_length=128, default='', unique=False)

    def __str__(self):
        return f'{self.name}-{self.version}'

    class Meta:
        unique_together = ('name', 'version')


class PDisDataset(models.Model):
    name = models.CharField(max_length=128, default='', unique=False)
    link = models.CharField(max_length=128, default='', unique=False)
    version = models.CharField(max_length=128, default='', unique=False)

    def __str__(self):
        return f'{self.name}-{self.version}'

    class Meta:
        unique_together = ('name', 'version')


class DrDiDataset(models.Model):
    name = models.CharField(max_length=128, default='', unique=False)
    link = models.CharField(max_length=128, default='', unique=False)
    version = models.CharField(max_length=128, default='', unique=False)

    def __str__(self):
        return f'{self.name}-{self.version}'

    class Meta:
        unique_together = ('name', 'version')


class ExpressionLevel(models.Model):
    tissue = models.ForeignKey('Tissue', on_delete=models.CASCADE)
    protein = models.ForeignKey('Protein', on_delete=models.CASCADE)
    expression_level = models.FloatField()

    class Meta:
        unique_together = ('tissue', 'protein')


class EnsemblGene(models.Model):
    name = models.CharField(max_length=15, unique=True)  # starts with ENSG...
    protein = models.ForeignKey('Protein', on_delete=models.CASCADE, related_name='ensg')


class Protein(models.Model):
    # According to https://www.uniprot.org/help/accession_numbers UniProt accession codes
    # are either 6 or 10 characters long

    uniprot_code = models.CharField(max_length=10)
    gene = models.CharField(max_length=128, default='')   # symbol
    protein_name = models.CharField(max_length=128, default='')
    entrez = models.CharField(max_length=128, default='')
    drugs = models.ManyToManyField('Drug', through='ProteinDrugInteraction',
                                   related_name='interacting_drugs')
    tissue_expression = models.ManyToManyField('Tissue', through='ExpressionLevel',
                                               related_name='interacting_drugs')

    class Meta:
        unique_together = ('uniprot_code', 'gene', 'entrez')

    def __str__(self):
        return self.gene


class Disorder(models.Model):
    mondo_id = models.CharField(max_length=7)
    label = models.CharField(max_length=256, default='')   # symbol
    icd10 = models.CharField(max_length=128, default='')
    proteins = models.ManyToManyField(
        'Protein', through='ProteinDisorderAssociation', related_name='associated_proteins')

    class Meta:
        unique_together = ('mondo_id', 'label', 'icd10')

    def __str__(self):
        return self.label


class ProteinDisorderAssociation(models.Model):
    pdis_dataset = models.ForeignKey(
        'PDisDataset', null=True, on_delete=models.CASCADE, related_name='pdis_dataset_relation')
    protein = models.ForeignKey('Protein', on_delete=models.CASCADE)
    disorder = models.ForeignKey('Disorder', on_delete=models.CASCADE)
    score = models.FloatField()

    class Meta:
        unique_together = ('pdis_dataset', 'protein', 'disorder')

    def __str__(self):
        return f'{self.pdis_dataset}-{self.protein}-{self.disorder}'


class DrugDisorderIndication(models.Model):
    drdi_dataset = models.ForeignKey(
        'DrDiDataset', null=True, on_delete=models.CASCADE, related_name='drdi_dataset_relation')
    drug = models.ForeignKey('Drug', on_delete=models.CASCADE)
    disorder = models.ForeignKey('Disorder', on_delete=models.CASCADE)

    class Meta:
        unique_together = ('drdi_dataset', 'drug', 'disorder')

    def __str__(self):
        return f'{self.drdi_dataset}-{self.drug}-{self.disorder}'


class Drug(models.Model):
    drug_id = models.CharField(max_length=10, unique=True)
    name = models.CharField(max_length=256, default='')
    status = models.CharField(max_length=128, default='')
    # in_trial = models.BooleanField(default=False)
    # in_literature = models.BooleanField(default=False)
    links = models.CharField(max_length=16*1024, default='')

    def __str__(self):
        return self.drug_id


class ProteinProteinInteraction(models.Model):
    ppi_dataset = models.ForeignKey(
        'PPIDataset', null=True, on_delete=models.CASCADE, related_name='ppi_dataset_relation')
    from_protein = models.ForeignKey('Protein', on_delete=models.CASCADE, related_name='interacting_proteins_out')
    to_protein = models.ForeignKey('Protein', on_delete=models.CASCADE, related_name='interacting_proteins_in')

    def validate_unique(self, exclude=None):
        p1p2_q = ProteinProteinInteraction.objects.filter(
            from_protein=self.from_protein,
            to_protein=self.to_protein,
            ppi_dataset=self.ppi_dataset
            )
        p2p1_q = ProteinProteinInteraction.objects.filter(
            from_protein=self.to_protein,
            to_protein=self.from_protein,
            ppi_dataset=self.ppi_dataset
            )

        if p1p2_q.exists() or p2p1_q.exists():
            raise ValidationError('Protein-Protein interaction must be unique!')

    def save(self, *args, **kwargs):
        self.validate_unique()
        super(ProteinProteinInteraction, self).save(*args, **kwargs)

    def __str__(self):
        return f'{self.ppi_dataset}-{self.from_protein}-{self.to_protein}'


class ProteinDrugInteraction(models.Model):
    pdi_dataset = models.ForeignKey(
        'PDIDataset', null=True, on_delete=models.CASCADE, related_name='pdi_dataset_relation')
    protein = models.ForeignKey('Protein', on_delete=models.CASCADE)
    drug = models.ForeignKey('Drug', on_delete=models.CASCADE)

    class Meta:
        unique_together = ('pdi_dataset', 'protein', 'drug')

    def __str__(self):
        return f'{self.pdi_dataset}-{self.protein}-{self.drug}'


class Task(models.Model):
    token = models.CharField(max_length=32, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    target = models.CharField(max_length=32, choices=[('drug', 'Drug'), ('drug-target', 'Drug Target')])

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
    nodes = models.TextField(null=True, default='')
    edges = models.TextField(null=True, default='')
    config = models.TextField(null=True, default='')
    groups = models.TextField(null=True, default='')
