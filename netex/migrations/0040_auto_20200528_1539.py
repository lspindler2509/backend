# Generated by Django 3.0.5 on 2020-05-28 13:39

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('netex', '0039_auto_20200528_1503'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='expressionlevel',
            unique_together={('tissue', 'protein')},
        ),
    ]
