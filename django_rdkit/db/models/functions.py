from __future__ import unicode_literals

import sys

from django.db import models

from django_rdkit.db.models.descriptors import DESCRIPTOR_MIXINS
from django_rdkit.db.models.fields import *

__all__ = []

module = sys.modules[__name__]

for mixin in DESCRIPTOR_MIXINS:
    F = type(str(mixin.descriptor_name.upper()), (mixin, models.Func,), {})
    setattr(module, F.__name__, F)
    __all__.append(F.__name__)


for name, function, fieldkls in [
        ('MOL', 'mol', MolField),
        ('MOL_FROM_SMILES', 'mol_from_smiles', MolField),
        ('MOL_FROM_SMARTS', 'mol_from_smarts', MolField),
        ('MOL_FROM_CTAB', 'mol_from_ctab', MolField),
        ('QMOL', 'qmol',  MolField),
        ('SMILES', 'mol_to_smiles', models.CharField),
        ('SMARTS', 'mol_to_smarts', models.CharField),
        ('CTAB', 'mol_to_ctab', models.TextField),
        ('INCHI', 'mol_inchi', models.TextField),
        ('INCHIKEY', 'mol_inchikey', models.TextField),
        ('FORMULA', 'mol_formula', models.TextField),
        ]:
    F = type(str(name), (models.Func,), 
             { 'function': function, 'output_field': fieldkls(),})
    setattr(module, F.__name__, F)
    __all__.append(F.__name__)


for validator in ['is_valid_smiles', 'is_valid_smarts', 'is_valid_ctab']:
    F = type(str(validator.upper()), (models.Func,), 
             { 'function': validator, 'output_field': models.BooleanField(),})
    setattr(module, F.__name__, F)
    __all__.append(F.__name__)


for fingerprint, fieldkls in [('morgan_fp', SfpField),
                              ('morganbv_fp', BfpField),
                              ('featmorgan_fp', SfpField),
                              ('featmorganbv_fp', BfpField),
                              ('rdkit_fp', BfpField),
                              ('atompair_fp', SfpField),
                              ('atompairbv_fp', BfpField),
                              ('torsion_fp', SfpField),
                              ('torsionbv_fp', BfpField),
                              ('layered_fp', BfpField),
                              ('maccs_fp', BfpField),
                              ('tanimoto_sml', models.FloatField),
                              ('dice_sml', models.FloatField),
                              ('tversky_sml', models.FloatField),
                              ('tanimoto_dist', models.FloatField),
                              ('dice_dist', models.FloatField),
                          ]:
    F = type(str(fingerprint.upper()), (models.Func,), 
             { 'function': fingerprint, 'output_field': fieldkls(),})
    setattr(module, F.__name__, F)
    __all__.append(F.__name__)


