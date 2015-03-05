from __future__ import unicode_literals

from django import VERSION as _django_version

from django.db.models import *

# Chem aggregate functions
# future -> from django_rdkit.models.aggregates import * 

from django_rdkit.models.fields import *

if _django_version[:2] > (1, 7):
    from django_rdkit.models.functions import *
