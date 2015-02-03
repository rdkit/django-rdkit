import sys

from django.db import models

from django_rdkit.db.models.descriptors import DESCRIPTOR_MIXINS
from django_rdkit.db.models.fields import *

__all__ = [ 
    "MOL", "QMOL", 
    "MORGAN_FP", "MORGANBV_FP",
]

module = sys.modules[__name__]

for mixin in DESCRIPTOR_MIXINS:
    F = type(mixin.descriptor_name.upper(), (mixin, models.Func,), {})
    setattr(module, F.__name__, F)
    __all__.append(F.__name__)


class MOL(models.Func):
    # the default template, resolving to mol(%(expressions)s)
    # should also work since the type name can be used as a function
    # template = '%(expressions)s::mol'
    function = 'mol'
    output_field = MoleculeField()


class QMOL(models.Func):
    function = 'qmol'


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
                          ]:
    F = type(fingerprint.upper(), (models.Func,), 
             { 'function': fingerprint, 'output_field': fieldkls(),})
    setattr(module, F.__name__, F)
    __all__.append(F.__name__)


