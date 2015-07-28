from __future__ import unicode_literals

import sys

from django.db import models
from django.db.models.expressions import Expression

from django_rdkit.models.fields import *
from django_rdkit.models.fields import MOL_DESCRIPTOR_MIXINS
from django_rdkit.models.fields import RXN_DESCRIPTOR_MIXINS

__all__ = []

module = sys.modules[__name__]

class _Func(models.Func):

    def __init__(self, *args, **kwargs):
        if not 'output_field' in kwargs:
            kwargs['output_field'] = self.default_output_field
        super(_Func, self).__init__(*args, **kwargs)


for mixin in MOL_DESCRIPTOR_MIXINS:
    _F = type(str(mixin.descriptor_name.upper()), (mixin, _Func,), {})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


for function, fieldkls in [('mol', MolField),
                           #('mol_from_smiles', MolField),
                           #('mol_from_smarts', MolField),
                           #('mol_from_ctab', MolField),
                           ('qmol',  MolField),
                           #('qmol_from_smiles', MolField),
                           #('qmol_from_ctab', MolField),
                           ('mol_to_smiles', models.CharField),
                           ('mol_to_smarts', models.CharField),
                           ('mol_to_ctab', models.TextField),
                           ('mol_inchi', models.TextField),
                           ('mol_inchikey', models.TextField),
                           ('mol_formula', models.TextField),
                       ]:
    _F = type(str(function.upper()), (_Func,),
             { 'function': function, 'default_output_field': fieldkls(),})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


class ConstructorFunc(_Func):
    #template = '%(function)s(cstring(%(expressions)s))'

    def as_sql(self, compiler, connection, function=None, template=None):
        connection.ops.check_expression_support(self)
        sql_parts = []
        params = []
        for arg in self.source_expressions:
            arg_sql, arg_params = compiler.compile(arg)
            sql_parts.append(arg_sql)
            params.extend(arg_params)
        if function is None:
            self.extra['function'] = self.extra.get('function', self.function)
        else:
            self.extra['function'] = function
        sql_parts[0] = 'cstring(%s)' % sql_parts[0]
        self.extra['expressions'] = self.extra['field'] = self.arg_joiner.join(sql_parts)
        template = template or self.extra.get('template', self.template)
        return template % self.extra, params


for constructor in ['mol_from_smiles', 'mol_from_smarts', 'mol_from_ctab',
                    'qmol_from_smiles', 'qmol_from_ctab',
                ]:
    _F = type(str(constructor.upper()), (ConstructorFunc,),
              { 'function': constructor,
                'default_output_field': MolField(),})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)
    

class ValidatorFunc(_Func):
    template = '%(function)s(cstring(%(expressions)s))'


for validator in ['is_valid_smiles', 'is_valid_smarts', 'is_valid_ctab']:
    _F = type(str(validator.upper()), (ValidatorFunc,),
             { 'function': validator,
               'default_output_field': models.BooleanField(),})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


for mixin in RXN_DESCRIPTOR_MIXINS:
    _F = type(str(mixin.descriptor_name.upper()), (mixin, _Func,), {})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


for function, fieldkls in [('reaction', RxnField),
                           ('reaction_from_smiles', RxnField),
                           ('reaction_from_smarts', RxnField),
                           ('reaction_from_ctab', RxnField),
                           ('reaction_to_smiles', models.CharField),
                           ('reaction_to_smarts', models.CharField),
                           ('reaction_to_ctab', models.TextField),
                           ('reaction_difference_fp', SfpField),
                           ('reaction_structural_bfp', BfpField),
                       ]:
    _F = type(str(function.upper()), (_Func,),
             { 'function': function, 'default_output_field': fieldkls(),})
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
    _F = type(str(fingerprint.upper()), (_Func,),
             { 'function': fingerprint, 'default_output_field': fieldkls(),})
    setattr(module, _F.__name__, _F)
    __all__.append(_F.__name__)


class DistanceExpression(Expression):

    def __init__(self, lhs, connector, rhs):
        super(DistanceExpression, self).__init__(output_field=models.FloatField())
        self.lhs, self.rhs = [
            arg if hasattr(arg, 'resolve_expression') else models.F(arg)
            for arg in (lhs, rhs)
        ]
        self.connector = connector

    def get_source_expressions(self):
        return [self.lhs, self.rhs]

    def set_source_expressions(self, exprs):
        self.lhs, self.rhs = exprs

    def resolve_expression(self, query=None, 
                           allow_joins=True, reuse=None, summarize=False):
        c = self.copy()
        c.is_summary = summarize
        c.lhs = self.lhs.resolve_expression(query, allow_joins, reuse, summarize)
        c.rhs = self.rhs.resolve_expression(query, allow_joins, reuse, summarize)
        return c

    def as_sql(self, compiler, connection):
        lhs_sql, lhs_params = compiler.compile(self.lhs)
        rhs_sql, rhs_params = compiler.compile(self.rhs)
        expression_wrapper = '(%s)'
        sql = "%s %s %s" % (lhs_sql, self.connector, rhs_sql)
        return expression_wrapper % sql, lhs_params + rhs_params


class TANIMOTO_DIST(DistanceExpression):
    def __init__(self, lhs, rhs):
        super(TANIMOTO_DIST, self).__init__(lhs, '<%%>', rhs)


__all__.append('TANIMOTO_DIST')


class DICE_DIST(DistanceExpression):
    def __init__(self, lhs, rhs):
        super(DICE_DIST, self).__init__(lhs, '<#>', rhs)


__all__.append('DICE_DIST')


class FMCS(models.Aggregate):
    function = 'fmcs'
    def __init__(self, expression):
        super(FMCS, self).__init__(expression, output_field=models.CharField())


__all__.append('FMCS')

