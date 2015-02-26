from __future__ import unicode_literals

import sys

from django.db import models
from django.db.models.expressions import Expression

from django_rdkit.models.fields import *
from django_rdkit.models.fields import DESCRIPTOR_MIXINS

__all__ = []

module = sys.modules[__name__]

for mixin in DESCRIPTOR_MIXINS:
    _F = type(str(mixin.descriptor_name.upper()), (mixin, models.Func,), {})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


for function, fieldkls in [('mol', MolField),
                           ('mol_from_smiles', MolField),
                           ('mol_from_smarts', MolField),
                           ('mol_from_ctab', MolField),
                           ('qmol',  MolField),
                           ('mol_to_smiles', models.CharField),
                           ('mol_to_smarts', models.CharField),
                           ('mol_to_ctab', models.TextField),
                           ('mol_inchi', models.TextField),
                           ('mol_inchikey', models.TextField),
                           ('mol_formula', models.TextField),
                       ]:
    _F = type(str(function.upper()), (models.Func,), 
             { 'function': function, 'output_field': fieldkls(),})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


for validator in ['is_valid_smiles', 'is_valid_smarts', 'is_valid_ctab']:
    _F = type(str(validator.upper()), (models.Func,), 
             { 'function': validator, 'output_field': models.BooleanField(),})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


for function, fieldkls in [('reaction', RxnField),
                           ('reaction_from_smiles', RxnField),
                           ('reaction_from_smarts', RxnField),
                           ('reaction_from_ctab', RxnField),
                           ('reaction_to_smiles', models.CharField),
                           ('reaction_to_smarts', models.CharField),
                           ('reaction_to_ctab', models.TextField),
                       ]:
    _F = type(str(function.upper()), (models.Func,), 
             { 'function': function, 'output_field': fieldkls(),})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


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
                              #('tanimoto_dist', models.FloatField),
                              #('dice_dist', models.FloatField),
                          ]:
    _F = type(str(fingerprint.upper()), (models.Func,), 
             { 'function': fingerprint, 'output_field': fieldkls(),})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


class TANIMOTO_DIST(Expression):

    def __init__(self, lhs, rhs):
        lhs, rhs = [
            arg if hasattr(arg, 'resolve_expression') else models.F(arg)
            for arg in (lhs, rhs)
        ]
        super(TANIMOTO_DIST, self).__init__(lhs, '<%%>', rhs, 
                                            models.FloatField())


__all__.append('TANIMOTO_DIST')


class DICE_DIST(Expression):

    def __init__(self, lhs, rhs):
        lhs, rhs = [
            arg if hasattr(arg, 'resolve_expression') else models.F(arg)
            for arg in (lhs, rhs)
        ]
        super(DICE_DIST, self).__init__(lhs, '<#>', rhs, 
                                            models.FloatField())


__all__.append('DICE_DIST')
