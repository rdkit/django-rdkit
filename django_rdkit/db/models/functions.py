from django.db import models

from django_rdkit.db.models.fields import MoleculeField

__all__ = [ "MOL", "QMOL", ]

class MOL(models.Func):
    # the default template, resolving to mol(%(expressions)s)
    # should also work since the type name can be used as a function
    # template = '%(expressions)s::mol'
    function = 'mol'
    output_field = MoleculeField

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


