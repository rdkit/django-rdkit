# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django_rdkit.db.models.fields


class Migration(migrations.Migration):

    dependencies = [
        ('app', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='compound',
            name='ffp2',
            field=django_rdkit.db.models.fields.BfpField(null=True),
        ),
        migrations.AddField(
            model_name='compound',
            name='mfp2',
            field=django_rdkit.db.models.fields.BfpField(null=True),
        ),
        migrations.AddField(
            model_name='compound',
            name='torsionbv',
            field=django_rdkit.db.models.fields.BfpField(null=True),
        ),
    ]
