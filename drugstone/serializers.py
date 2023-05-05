# Serializers define the API representation.
import json
from rest_framework import serializers
from drugstone import models
from drugstone.models import Protein, Task, Drug, ProteinDrugInteraction, \
    Tissue, ProteinProteinInteraction, Network, ProteinDisorderAssociation, Disorder, DrugDisorderIndication


class PDIDatasetSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.PDIDataset
        fields = '__all__'


class PPIDatasetSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.PPIDataset
        fields = '__all__'


class PDisDatasetSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.PDisDataset
        fields = "__all__"


class DrDisDatasetSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.DrDiDataset
        fields = "__all__"


class ProteinNodeSerializer(serializers.ModelSerializer):
    drugstone_id = serializers.SerializerMethodField()
    uniprot = serializers.SerializerMethodField()
    symbol = serializers.SerializerMethodField()
    ensg = serializers.SerializerMethodField()
    entrez = serializers.SerializerMethodField()

    def get_drugstone_id(self, obj):
        return [f"p{obj.id}"]

    def get_uniprot(self, obj):
        return [obj.uniprot_code]

    def get_symbol(self, obj):
        return [obj.gene]

    def get_entrez(self, obj):
        return [obj.entrez]

    def get_ensg(self, obj) -> str:
        """Since ENSG has a many to one relationship to the Protein table,
        return a list of all matching ensg names.

        Args:
            obj (Protein): Protein object

        Returns:
            str: list of all matching ENSG numbers
        """
        return [x.name for x in obj.ensg.all()]

    class Meta:
        model = Protein
        fields = ["drugstone_id", "uniprot", "symbol", "protein_name", "entrez", "ensg"]


class ProteinSerializer(serializers.ModelSerializer):
    drugstone_id = serializers.SerializerMethodField()
    uniprot = serializers.SerializerMethodField()
    symbol = serializers.SerializerMethodField()
    ensg = serializers.SerializerMethodField()

    def get_drugstone_id(self, obj):
        return f"p{obj.id}"

    def get_uniprot(self, obj):
        return obj.uniprot_code

    def get_symbol(self, obj):
        return obj.gene

    def get_ensg(self, obj) -> str:
        """Since ENSG has a many to one relationship to the Protein table,
        return a list of all matching ensg names.

        Args:
            obj (Protein): Protein object

        Returns:
            str: list of all matching ENSG numbers
        """
        return [x.name for x in obj.ensg.all()]

    class Meta:
        model = Protein
        fields = ["drugstone_id", "uniprot", "symbol", "protein_name", "entrez", "ensg"]


class DrugSerializer(serializers.ModelSerializer):
    drugstone_id = serializers.SerializerMethodField()
    trial_links = serializers.SerializerMethodField()
    label = serializers.SerializerMethodField()

    def get_drugstone_id(self, obj):
        return f"dr{obj.id}"

    def get_trial_links(self, obj):
        return [] if obj.links == "" else obj.links.split(";")

    def get_label(self, obj):
        return obj.name

    class Meta:
        model = Drug
        fields = ["drugstone_id", "drug_id", "label", "status", "trial_links"]


class DisorderSerializer(serializers.ModelSerializer):
    drugstone_id = serializers.SerializerMethodField()
    icd_10 = serializers.SerializerMethodField()
    disorder_id = serializers.SerializerMethodField()

    def get_drugstone_id(self, obj):
        return f"di{obj.id}"

    def get_icd_10(self, obj):
        return obj.icd10[1: len(obj.icd10) - 1].split(",")

    def get_disorder_id(self, obj):
        return obj.mondo_id

    class Meta:
        model = Disorder
        fields = ["drugstone_id", "label", "icd_10", "disorder_id"]


class ProteinProteinInteractionSerializer(serializers.ModelSerializer):
    dataset = serializers.SerializerMethodField()
    protein_a = serializers.SerializerMethodField()
    protein_b = serializers.SerializerMethodField()

    def get_dataset(self, obj):
        return obj.ppi_dataset.name

    def get_protein_a(self, obj):
        return f"p{obj.from_protein.id}"

    def get_protein_b(self, obj):
        return f"p{obj.to_protein.id}"

    class Meta:
        model = ProteinProteinInteraction
        fields = ["dataset", "protein_a", "protein_b"]


class ProteinDrugInteractionSerializer(serializers.ModelSerializer):
    dataset = serializers.SerializerMethodField()
    protein = serializers.SerializerMethodField()
    drug = serializers.SerializerMethodField()
    actions = serializers.SerializerMethodField()

    def get_dataset(self, obj):
        return obj.pdi_dataset.name

    def get_protein(self, obj):
        return f"p{obj.protein.id}"

    def get_drug(self, obj):
        return f"dr{obj.drug.id}"

    def get_actions(self, obj):
        if obj.actions:
            return json.loads(obj.actions)
        return []

    class Meta:
        model = ProteinDrugInteraction
        fields = ['dataset', 'protein', 'drug', 'actions']


class ProteinDisorderAssociationSerializer(serializers.ModelSerializer):
    dataset = serializers.SerializerMethodField()
    protein = serializers.SerializerMethodField()
    disorder = serializers.SerializerMethodField()
    score = serializers.SerializerMethodField()

    def get_dataset(self, obj):
        return obj.pdis_dataset.name

    def get_protein(self, obj):
        return f"p{obj.protein.id}"

    def get_disorder(self, obj):
        return f"di{obj.disorder.id}"

    def get_score(self, obj):
        return float(obj.score)

    class Meta:
        model = ProteinDisorderAssociation
        fields = ["dataset", "protein", "disorder", "score"]


class DrugDisorderIndicationSerializer(serializers.ModelSerializer):
    dataset = serializers.SerializerMethodField()
    drug = serializers.SerializerMethodField()
    disorder = serializers.SerializerMethodField()

    def get_dataset(self, obj):
        return obj.drdi_dataset.name

    def get_drug(self, obj):
        return f"dr{obj.drug.id}"

    def get_disorder(self, obj):
        return f"di{obj.disorder.id}"

    class Meta:
        model = DrugDisorderIndication
        fields = ["dataset", "drug", "disorder"]


class TaskSerializer(serializers.ModelSerializer):
    parameters = serializers.SerializerMethodField()

    def get_parameters(self, obj):
        return json.loads(obj.parameters)

    class Meta:
        model = Task
        fields = [
            "algorithm",
            "target",
            "parameters",
            "job_id",
            "worker_id",
            "progress",
            "status",
            "created_at",
            "started_at",
            "finished_at",
            "done",
            "failed",
        ]


class NetworkSerializer(serializers.ModelSerializer):
    # nodes = serializers.SerializerMethodField()
    # edges = serializers.SerializerMethodField()
    # config = serializers.SerializerMethodField()

    class Meta:
        model = Network
        fields = "__all__"


#    def get_nodes(self,obj):
#        return json.loads(obj.nodes)
#
#    def get_edges(self,obj):
#        return json.loads(obj.edges)
#
#    def get_config(self,obj):
#        return json.loads(obj.config)


class TaskStatusSerializer(serializers.ModelSerializer):
    class Meta:
        model = Task
        fields = [
            "algorithm",
            "target",
            "progress",
            "status",
            "created_at",
            "started_at",
            "finished_at",
            "done",
            "failed",
        ]


class TissueSerializer(serializers.ModelSerializer):
    drugstone_id = serializers.SerializerMethodField()

    def get_drugstone_id(self, obj):
        return f"{obj.id}"

    class Meta:
        model = Tissue
        fields = ["drugstone_id", "name"]
