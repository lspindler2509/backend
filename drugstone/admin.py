from django.contrib import admin
import drugstone.models as models

admin.site.register(models.Protein)
admin.site.register(models.Tissue)
admin.site.register(models.Drug)
admin.site.register(models.Task)

admin.site.register(models.PPIDataset)
admin.site.register(models.ProteinProteinInteraction)

admin.site.register(models.PDIDataset)
admin.site.register(models.ProteinDrugInteraction)
