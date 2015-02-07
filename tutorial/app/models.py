from django_rdkit.db import models


class Compound(models.Model):

    name = models.CharField(max_length=256)
    molecule = models.MolField()

    torsionbv = models.BfpField(null=True)
    mfp2 = models.BfpField(null=True)
    ffp2 = models.BfpField(null=True)

