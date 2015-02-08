# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django_rdkit.db.models.fields


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Compound',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=256)),
                ('molecule', django_rdkit.db.models.fields.MolField()),
                ('torsionbv', django_rdkit.db.models.fields.BfpField(null=True)),
                ('mfp2', django_rdkit.db.models.fields.BfpField(null=True)),
                ('ffp2', django_rdkit.db.models.fields.BfpField(null=True)),
            ],
        ),
    ]
