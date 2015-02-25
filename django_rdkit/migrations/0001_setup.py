# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations
from django_rdkit.operations import RDKitExtension

class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        RDKitExtension(),
    ]
