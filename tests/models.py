#from __future__ import unicode_literals

from django_rdkit.db import models

class TestModel(models.Model):

    molecule = models.MoleculeField()
