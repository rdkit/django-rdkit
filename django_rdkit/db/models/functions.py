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

    def __init__(self, *expressions, **extra):
      if len(expressions) != 1:
          raise ValueError('expressions must have exacly 1 element')
      super(MOL, self).__init__(*expressions, **extra)


class QMOL(models.Func):
    function = 'qmol'

    def __init__(self, *expressions, **extra):
      if len(expressions) != 1:
          raise ValueError('expressions must have exacly 1 element')
      super(QMOL, self).__init__(*expressions, **extra)


class MORGAN_FP(models.Func):
    function = 'morgan_fp'
    output_field = SfpField()

    def __init__(self, *expressions, **extra):
      if len(expressions) > 2:
          raise ValueError('expressions must have at most 2 elements')
      super(MORGAN_FP, self).__init__(*expressions, **extra)


class MORGANBV_FP(models.Func):
    function = 'morganbv_fp'
    output_field = BfpField()

    def __init__(self, *expressions, **extra):
      if len(expressions) > 2:
          raise ValueError('expressions must have at most 2 elements')
      super(MORGANBV_FP, self).__init__(*expressions, **extra)


