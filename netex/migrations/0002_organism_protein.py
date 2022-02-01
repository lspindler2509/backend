# Generated by Django 3.0.4 on 2020-03-28 13:41

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('netex', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Organism',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=128)),
            ],
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('uniprot_code', models.CharField(max_length=10)),
                ('protein_name', models.CharField(max_length=128)),
                ('gene_name', models.CharField(max_length=128)),
                ('organism', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='netex.Organism')),
                ('protein_group', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='netex.ProteinGroup')),
            ],
        ),
    ]
