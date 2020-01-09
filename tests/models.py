from django_rdkit import models

class MoleculeModel(models.Model):

    molecule = models.MolField(null=True)


class ReactionModel(models.Model):

    rxn = models.RxnField()


class SfpModel(models.Model):

    sfp = models.SfpField(null=True)


class BfpModel(models.Model):

    bfp = models.BfpField(null=True)


class SmilesModel(models.Model):

    smiles = models.CharField(max_length=2048, blank=True, null=False)
    molecule = models.MolField(null=True)


class SmartsModel(models.Model):

    smarts = models.CharField(max_length=2048, blank=True, null=False)


class CtabModel(models.Model):

    ctab = models.TextField(blank=True, null=False)


#class TestModel(models.Model):#
#
#    molecule = models.MolField()
#    sfp = models.SfpField()
#    bfb = models.BfpField()
