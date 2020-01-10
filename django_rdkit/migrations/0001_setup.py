# -*- coding: utf-8 -*-
from django.db import migrations
from django_rdkit.operations import RDKitExtension

class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        RDKitExtension(),
    ]
