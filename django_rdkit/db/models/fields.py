from __future__ import unicode_literals

from django.utils import six
from django.utils.translation import ugettext_lazy as _
from django.db.models import Lookup, Transform
from django.db.models.fields import *

from rdkit import Chem
from rdkit.Chem import Mol
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect, SparseIntVect


from django_rdkit.db.models.descriptors import DESCRIPTOR_MIXINS


__all__ = ["MolField", "BfpField", "SfpField",]
 

class ChemField(Field):

    def __init__(self, verbose_name=None, chem_index=True, *args, **kwargs):
        self.chem_index = chem_index
        kwargs['verbose_name'] = verbose_name
        super(ChemField, self).__init__(*args, **kwargs)
    
    #def select_format(self, compiler, sql, params):
    #    """
    #    Custom format for select clauses. For example, GIS columns need to be
    #    selected as AsText(table.col) on MySQL as the table.col data 
    #    can't be used by Django.
    #    """
    #    return sql, params

    def deconstruct(self):
        name, path, args, kwargs = super(ChemField, self).deconstruct()
        # include chem_index if not the default value.
        if self.chem_index is not True:
            kwargs['chem_index'] = self.chem_index
        return name, path, args, kwargs

    #def get_placeholder(self, value, compiler, connection):
    #    """
    #    Returns the placeholder for the geometry column for the
    #    given value.
    #    """
    #    if hasattr(value, 'as_sql'):
    #        # No geometry value used for F expression, substitute in
    #        # the column name instead.
    #        sql, _ = compiler.compile(value)
    #        return sql
    #    else:
    #        return '%s'


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

    def from_db_value(self, value, connection, context):
        if value is None:
            return value
        return Mol(bytes(value))

    def to_python(self, value):
        if value is None or isinstance(value, Mol):
            return value
        elif isinstance(value, six.string_types):
            # The string case. A SMILES is assumed.
            return Chem.MolFromSmiles(str(value))
        elif isinstance(value, six.buffer_types):
            return Mol(bytes(value))
        else:
            raise ValidationError("Invalid input for a Mol instance")

    def get_prep_value(self, value):
        # convert the Molecule instance to the value used by the 
        # db driver
        if isinstance(value, six.string_types):
            # The string case. A SMILES is assumed.
            value = Chem.MolFromSmiles(str(value))
        if isinstance(value, Mol):
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

class DescriptorTransform(Transform):

    def as_sql(self, qn, connection):
        lhs, params = qn.compile(self.lhs)
        return "%s(%s)" % (self.function, lhs), params
    

DESCRIPTOR_TRANFORMS = [
    type(str('{0}_Transform'.format(mixin.descriptor_name.upper())),
         (mixin, DescriptorTransform,),
         { 'lookup_name': mixin.descriptor_name.upper(), }
     )
    for mixin in DESCRIPTOR_MIXINS
]


for Transform in DESCRIPTOR_TRANFORMS:
    MolField.register_lookup(Transform)


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

    def from_db_value(self, value, connection, context):
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






