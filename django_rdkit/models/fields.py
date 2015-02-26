from __future__ import unicode_literals

from django.utils import six
from django.utils.translation import ugettext_lazy as _
from django.db.models import Lookup, Transform
from django.db.models.fields import *

from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect, SparseIntVect


__all__ = ["MolField", "RxnField", "BfpField", "SfpField",]
 

class ChemField(Field):

    def __init__(self, verbose_name=None, chem_index=True, *args, **kwargs):
        self.chem_index = chem_index
        kwargs['verbose_name'] = verbose_name
        super(ChemField, self).__init__(*args, **kwargs)
    
    def deconstruct(self):
        name, path, args, kwargs = super(ChemField, self).deconstruct()
        # include chem_index if not the default value.
        if self.chem_index is not True:
            kwargs['chem_index'] = self.chem_index
        return name, path, args, kwargs


##########################################
# Molecule Field

class MolField(ChemField):

    description = _("Molecule")

    def db_type(self, connection):
        return 'mol'
    
    def get_placeholder(self, value, compiler, connection):
        if hasattr(value, 'as_sql'):
            # No value used for expressions, substitute in
            # the column name instead.
            sql, _ = compiler.compile(value)
            return sql
        else:
            return 'mol_from_pkl(%s)'

    def select_format(self, compiler, sql, params):
        return 'mol_to_pkl(%s)' % sql, params

    def from_db_value(self, value, expression, connection, context):
        if value is None:
            return value
        return Chem.Mol(bytes(value))

    def to_python(self, value):
        if value is None or isinstance(value, Chem.Mol):
            return value
        elif isinstance(value, six.string_types):
            # The string case. A SMILES is assumed.
            return Chem.MolFromSmiles(str(value))
        elif isinstance(value, six.buffer_types):
            return Chem.Mol(bytes(value))
        else:
            raise ValidationError("Invalid input for a Mol instance")

    def get_prep_value(self, value):
        # convert the Molecule instance to the value used by the 
        # db driver
        if isinstance(value, six.string_types):
            # The string case. A SMILES is assumed.
            value = Chem.MolFromSmiles(str(value))
        if isinstance(value, Chem.Mol):
            value = six.memoryview(value.ToBinary())
        return value

    # don't reimplement db-specific preparation of query values for now
    # def get_db_prep_value(self, value, connection, prepared=False):
    #    return value

    def get_prep_lookup(self, lookup_type, value):
        "Perform preliminary non-db specific lookup checks and conversions"
        supported_lookup_types = (
            ['hassubstruct', 'issubstruct', 'exact',] +
            [T.lookup_name for T in DESCRIPTOR_TRANFORMS]
        )
        if lookup_type in supported_lookup_types:
            return value
        raise TypeError("Field has invalid lookup: %s" % lookup_type)

    # this will be probably needed.
    #def get_db_prep_lookup(lookup_type, value, connection, prepared=False):
    #    if not prepared:
    #        value = self.get_prep_lookup(lookup_type, value)
    #    return value


##########################################
# Reaction Field

class RxnField(ChemField):

    description = _("Reaction")

    def db_type(self, connection):
        return 'reaction'
    
    #def get_placeholder(self, value, compiler, connection):
    #    if hasattr(value, 'as_sql'):
    #        # No value used for expressions, substitute in
    #        # the column name instead.
    #        sql, _ = compiler.compile(value)
    #        return sql
    #    else:
    #        return 'reaction_from_pkl(%s)'

    #def select_format(self, compiler, sql, params):
    #    return 'reaction_to_pkl(%s)' % sql, params

    def from_db_value(self, value, expression, connection, context):
        if value is None:
            return value
        return Chem.ReactionFromSmarts(value, useSmiles=True)
        #return Chem.ChemicalReaction(bytes(value))

    def to_python(self, value):
        if value is None or isinstance(value, Chem.ChemicalReaction):
            return value
        elif isinstance(value, six.string_types):
            # The string case. A reaction SMILES is expected.
            return Chem.ReactionFromSmarts(str(value), useSmiles=True)
        #elif isinstance(value, six.buffer_types):
        #    return Chem.ChemicalReaction(bytes(value))
        else:
            raise ValidationError("Invalid input for a ChemicalReaction instance")

    def get_prep_value(self, value):
        # convert the ChemicalReaction instance to the value used by the 
        # db driver
        #if isinstance(value, six.string_types):
        #    # The string case. A reaction SMILES is assumed.
        #    value = Chem.ReactionFromSmarts(str(value), useSmiles=True)
        if isinstance(value, Chem.ChemicalReaction):
            #value = six.memoryview(value.ToBinary())
            value = Chem.ReactionToSmiles(value)
        return value

    # don't reimplement db-specific preparation of query values for now
    # def get_db_prep_value(self, value, connection, prepared=False):
    #    return value

    def get_prep_lookup(self, lookup_type, value):
        "Perform preliminary non-db specific lookup checks and conversions"
        supported_lookup_types = (
        )
        if lookup_type in supported_lookup_types:
            return value
        raise TypeError("Field has invalid lookup: %s" % lookup_type)

    # this will be probably needed.
    #def get_db_prep_lookup(lookup_type, value, connection, prepared=False):
    #    if not prepared:
    #        value = self.get_prep_lookup(lookup_type, value)
    #    return value


########################################################
# Binary Fingerprint Field

class BfpField(ChemField):

    description = _("Binary Fingerprint")

    def db_type(self, connection):
        return 'bfp'
    
    def get_placeholder(self, value, compiler, connection):
        if hasattr(value, 'as_sql'):
            # No value used for expressions, substitute in
            # the column name instead.
            sql, _ = compiler.compile(value)
            return sql
        else:
            return 'bfp_from_binary_text(%s)'

    def select_format(self, compiler, sql, params):
        return 'bfp_to_binary_text(%s)' % sql, params

    def from_db_value(self, value, expression, connection, context):
        if value is None:
            return value
        return DataStructs.CreateFromBinaryText(bytes(value))

    def to_python(self, value):
        if value is None or isinstance(value, ExplicitBitVect):
            return value
        elif isinstance(value, six.buffer_types):
            return DataStructs.CreateFromBinaryText(bytes(value))
        else:
            raise ValidationError("Invalid input for a Bfp instance")

    def get_prep_value(self, value):
        # convert the ExplicitBitVect instance to the value used by the 
        # db driver
        if isinstance(value, ExplicitBitVect):
            value = six.memoryview(DataStructs.BitVectToBinaryText(value))
        return value

    def get_prep_lookup(self, lookup_type, value):
        if lookup_type in [
                'lt', 'lte', 'exact', 'gte', 'gt', 'ne', 
                'tanimoto', 'dice']:
            return value
        raise TypeError("Field has invalid lookup: %s" % lookup_type)

    #def get_db_prep_lookup(lookup_type, value, connection, prepared=False):
    #    if not prepared:
    #        value = self.get_prep_lookup(lookup_type, value)
    #    return value


########################################################
# Sparse Integer Vector Fingerprint Field

class SfpField(ChemField):

    description = _("Sparse Integer Vector Fingerprint")

    def db_type(self, connection):
        return 'sfp'
    
    #def to_python(self, value):
    #    return value

    #def get_prep_value(self, value):
    #    return value

    def get_prep_lookup(self, lookup_type, value):
        if lookup_type in [
                'lt', 'lte', 'exact', 'gte', 'gt', 'ne', 
                'tanimoto', 'dice']:
            return value
        raise TypeError("Field has invalid lookup: %s" % lookup_type)

    #def get_db_prep_lookup(lookup_type, value, connection, prepared=False):
    #    if not prepared:
    #        value = self.get_prep_lookup(lookup_type, value)
    #    return value


###############################################################
# MolField lookup operations, substruct and exact searches

class HasSubstruct(Lookup):

    lookup_name = 'hassubstruct'

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s @> %s' % (lhs, rhs), params

MolField.register_lookup(HasSubstruct)


class IsSubstruct(Lookup):

    lookup_name = 'issubstruct'

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s <@ %s' % (lhs, rhs), params

MolField.register_lookup(IsSubstruct)


class SameStructure(Lookup):

    lookup_name = 'exact'

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        #return '%s @= %s' % (lhs, rhs), params
        return '%s <@ %s AND %s @> %s' % (lhs, rhs, lhs, rhs), params + params

MolField.register_lookup(SameStructure)


##########################################
# MolField transforms and descriptors

DESCRIPTORS = [
    ('hba', IntegerField),
    ('hbd', IntegerField),
    ('numatoms', IntegerField),
    ('numheavyatoms', IntegerField),
    ('numrotatablebonds', IntegerField),
    ('numheteroatoms', IntegerField),
    ('numrings', IntegerField),
    ('numaromaticrings', IntegerField),
    ('numaliphaticrings', IntegerField),
    ('numsaturatedrings', IntegerField),
    ('numaromaticheterocycles', IntegerField),
    ('numaliphaticheterocycles', IntegerField),
    ('numsaturatedheterocycles', IntegerField),
    ('numaromaticcarbocycles', IntegerField),
    ('numaliphaticcarbocycles', IntegerField),
    ('numsaturatedcarbocycles', IntegerField),
    ('amw', FloatField),
    ('logp', FloatField),
    ('tpsa', FloatField),
    ('fractioncsp3', FloatField),
    ('chi0v', FloatField),
    ('chi1v', FloatField),
    ('chi2v', FloatField),
    ('chi3v', FloatField),
    ('chi5v', FloatField),
    ('chi0n', FloatField),
    ('chi1n', FloatField),
    ('chi2n', FloatField),
    ('chi3n', FloatField),
    ('chi5n', FloatField),
    ('kappa1', FloatField),
    ('kappa2', FloatField),
    ('kappa3', FloatField),
    ('kappa4', FloatField),
    ('murckoscaffold', MolField),
]


def make_mixin(name, field):
    return type(
        str('{0}_Mixin'.format(name.upper())), 
        (object,),
        { 
            'descriptor_name': name, 
            'function': 'mol_{0}'.format(name),
            'output_field': field, 
        },
    )


DESCRIPTOR_MIXINS = [
    make_mixin(d, fieldkls()) for d, fieldkls in DESCRIPTORS 
]


class DescriptorTransform(Transform):

    def as_sql(self, qn, connection):
        lhs, params = qn.compile(self.lhs)
        return "%s(%s)" % (self.function, lhs), params
    

DESCRIPTOR_TRANFORMS = [
    type(str('{0}_Transform'.format(mixin.descriptor_name.upper())),
         (mixin, DescriptorTransform,),
         { 'lookup_name': mixin.descriptor_name, }
     )
    for mixin in DESCRIPTOR_MIXINS
]


for Transform in DESCRIPTOR_TRANFORMS:
    MolField.register_lookup(Transform)


####################################################################
# Fingerprint Fields lookup operations, similarity searches

class TanimotoSimilar(Lookup):

    lookup_name = 'tanimoto'

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s %%%% %s' % (lhs, rhs), params


BfpField.register_lookup(TanimotoSimilar)
SfpField.register_lookup(TanimotoSimilar)


class DiceSimilar(Lookup):

    lookup_name = 'dice'

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s # %s' % (lhs, rhs), params


BfpField.register_lookup(DiceSimilar)
SfpField.register_lookup(DiceSimilar)


class NotEqual(Lookup):

    lookup_name = 'ne'

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s <> %s' % (lhs, rhs), params


BfpField.register_lookup(NotEqual)
SfpField.register_lookup(NotEqual)






