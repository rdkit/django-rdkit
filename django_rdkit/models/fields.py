from enum import Enum

from django.utils.translation import gettext_lazy as _
from django.db.models import Lookup, Transform, Func, Value
from django.db.models.fields import *
from django.core.exceptions import ValidationError, ImproperlyConfigured
from django import forms
from django.conf import settings

from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect


__all__ = ["MolField", "RxnField", "BfpField", "SfpField",]

##########################################
# Molecule Field

class MolSerialization(Enum):
    BINARY = 'BINARY'
    TEXT = 'TEXT'

try:
    MOL_SERIALIZATION = MolSerialization(
        getattr(settings, 'DJANGO_RDKIT_MOL_SERIALIZATION', 'BINARY')
        )
except ValueError:
    raise ImproperlyConfigured(
        'An invalid DJANGO_RDKIT_MOL_SERIALIZATION value was found in the settings. '
        f'Supported values are {[v.value for v in MolSerialization]}'
        )


class MolFieldPklMixin:

    def get_placeholder(self, value, compiler, connection):
        # define if/how the value assigned to this field is
        # to be wrapped into an insertion query
        if hasattr(value, 'as_sql'):
            return '%s'
        else:
            return 'mol_from_pkl(%s)'

    def select_format(self, compiler, sql, params):
        # format to use when the corresponding column appears
        # in select clauses.
        return 'mol_to_pkl(%s)' % sql, params

    def from_db_value(self, value, expression, connection):
        # convert the value returned by the database driver
        # into the desired Python data type
        if value is None:
            return value
        return Chem.Mol(bytes(value))

    def get_prep_value(self, value):
        # convert from Python to the value to be used in queries
        if isinstance(value, str):
            value = self.text_to_mol(value)
        if isinstance(value, Chem.Mol):
            value = memoryview(value.ToBinary())
        return value


class MolFieldSmilesMixin:

    def from_db_value(self, value, expression, connection):
        # convert the value returned by the database driver
        # into the desired Python data type
        if value is None:
            return value
        return Chem.MolFromSmiles(value)

    def get_prep_value(self, value):
        # convert from Python to the value to be used in queries
        if isinstance(value, str):
            value = self.text_to_mol(value)
        if isinstance(value, Chem.Mol):
            value = Chem.MolToSmiles(value)
        return value


MolFieldSerializationMixin = {
    MolSerialization.BINARY: MolFieldPklMixin,
    MolSerialization.TEXT: MolFieldSmilesMixin,
}[MOL_SERIALIZATION]


class MolField(MolFieldSerializationMixin, Field):

    description = _("Molecule")

    def db_type(self, connection):
        # return the database column data type for this field
        return 'mol'

    def to_python(self, value):
        # convert the input value into the expected Python data
        # types (called during input cleanup and prior to field
        # validation)
        if value is None or isinstance(value, Chem.Mol):
            return value
        elif isinstance(value, str):
            return self.text_to_mol(value)
        elif isinstance(value, (bytes, bytearray, memoryview)):
            return Chem.Mol(bytes(value))
        else:
            raise ValidationError("Invalid input for a Mol instance")

    @staticmethod
    def text_to_mol(value):
        value = str(value)
        mol = (Chem.MolFromSmiles(value)
            or Chem.MolFromMolBlock(value)
            or Chem.inchi.MolFromInchi(value))
        if mol is None:
            raise ValidationError("Invalid input for a Mol instance")
        return mol

    def get_prep_lookup(self, lookup_type, value):
        "Perform preliminary non-db specific lookup checks and conversions"
        supported_lookup_types = (
            ['hassubstruct', 'issubstruct', 'exact', 'isnull',] +
            [T.lookup_name for T in MOL_DESCRIPTOR_TRANFORMS]
        )
        if lookup_type in supported_lookup_types:
            return value
        raise TypeError("Field has invalid lookup: %s" % lookup_type)

    def formfield(self, **kwargs):
        # Use TextField as default input form to accommodate line breaks needed for molBlocks
        defaults = {'form_class': forms.CharField, 'strip': False, 'widget':forms.Textarea}
        defaults.update(kwargs)
        return super().formfield(**defaults)

##########################################
# Reaction Field

class RxnField(Field):

    description = _("Reaction")

    def db_type(self, connection):
        return 'reaction'

    def from_db_value(self, value, expression, connection):
        if value is None:
            return value
        return Chem.ReactionFromSmarts(value, useSmiles=True)

    def to_python(self, value):
        if value is None or isinstance(value, Chem.ChemicalReaction):
            return value
        elif isinstance(value, str):
            # The string case. A reaction SMILES is expected.
            return Chem.ReactionFromSmarts(str(value), useSmiles=True)
        else:
            raise ValidationError("Invalid input for a ChemicalReaction instance")

    def get_prep_value(self, value):
        if isinstance(value, Chem.ChemicalReaction):
            value = Chem.ReactionToSmiles(value)
        return value

    def get_prep_lookup(self, lookup_type, value):
        "Perform preliminary non-db specific lookup checks and conversions"
        supported_lookup_types = (
            ['hassubstruct', 'issubstruct', 'isnull',] + #'exact',] +
            [T.lookup_name for T in RXN_DESCRIPTOR_TRANFORMS]
        )
        if lookup_type in supported_lookup_types:
            return value
        raise TypeError("Field has invalid lookup: %s" % lookup_type)


########################################################
# Binary Fingerprint Field

class BfpField(Field):

    description = _("Binary Fingerprint")

    def db_type(self, connection):
        return 'bfp'

    def get_placeholder(self, value, compiler, connection):
        if hasattr(value, 'as_sql'):
            return '%s'
        else:
            return 'bfp_from_binary_text(%s)'

    def select_format(self, compiler, sql, params):
        return 'bfp_to_binary_text(%s)' % sql, params

    def from_db_value(self, value, expression, connection):
        if value is None:
            return value
        return DataStructs.CreateFromBinaryText(bytes(value))

    def to_python(self, value):
        if value is None or isinstance(value, ExplicitBitVect):
            return value
        elif isinstance(value, (bytes, bytearray, memoryview)):
            return DataStructs.CreateFromBinaryText(bytes(value))
        else:
            raise ValidationError("Invalid input for a Bfp instance")

    def get_prep_value(self, value):
        # convert the ExplicitBitVect instance to the value used by the
        # db driver
        if isinstance(value, ExplicitBitVect):
            value = memoryview(DataStructs.BitVectToBinaryText(value))
        return value

    def get_prep_lookup(self, lookup_type, value):
        if lookup_type in [
                'lt', 'lte', 'exact', 'isnull', 'gte', 'gt', 'ne',
                'tanimoto', 'dice']:
            return value
        raise TypeError("Field has invalid lookup: %s" % lookup_type)


########################################################
# Sparse Integer Vector Fingerprint Field

class SfpField(Field):

    description = _("Sparse Integer Vector Fingerprint")

    def db_type(self, connection):
        return 'sfp'

    def get_prep_lookup(self, lookup_type, value):
        if lookup_type in [
                'lt', 'lte', 'exact', 'isnull', 'gte', 'gt', 'ne',
                'tanimoto', 'dice']:
            return value
        raise TypeError("Field has invalid lookup: %s" % lookup_type)


###################################################################
# MolField/RxnField lookup operations, substruct and exact searches

class MolLookupMixin:

    def get_prep_lookup(self):
        if self.rhs_is_direct_value():
            if isinstance(self.rhs, Chem.Mol):
                if MOL_SERIALIZATION == MolSerialization.BINARY:
                    self.rhs = self.rhs.ToBinary()
                    self.rhs = Func(self.rhs, function='mol_from_pkl')
                elif MOL_SERIALIZATION == MolSerialization.TEXT:
                    self.rhs = Value(Chem.MolToSmiles(self.rhs))
                else:
                    # this should never happen, because
                    # MOL_SERIALIZATION is validated at
                    # import time
                    raise NotImplementedError
            else:
                self.rhs = Value(self.rhs)
        return super().get_prep_lookup()


class HasSubstruct(Lookup):

    lookup_name = 'hassubstruct'
    prepare_rhs = True

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s @> %s' % (lhs, rhs), params


class HasMolSubstruct(MolLookupMixin, HasSubstruct):
    pass

MolField.register_lookup(HasMolSubstruct)
RxnField.register_lookup(HasSubstruct)


class HasSubstructFP(Lookup):

    lookup_name = 'hassubstructfp'

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s ?> %s' % (lhs, rhs), params


RxnField.register_lookup(HasSubstructFP)


class IsSubstruct(Lookup):

    lookup_name = 'issubstruct'
    prepare_rhs = True

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s <@ %s' % (lhs, rhs), params

class IsMolSubstruct(MolLookupMixin, IsSubstruct):
    pass

MolField.register_lookup(IsMolSubstruct)
RxnField.register_lookup(IsSubstruct)


class IsSubstructFP(Lookup):

    lookup_name = 'issubstructfp'

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        return '%s ?< %s' % (lhs, rhs), params

RxnField.register_lookup(IsSubstructFP)


class SameStructure(MolLookupMixin, Lookup):

    lookup_name = 'exact'
    prepare_rhs = True

    def as_sql(self, qn, connection):
        lhs, lhs_params = self.process_lhs(qn, connection)
        rhs, rhs_params = self.process_rhs(qn, connection)
        params = lhs_params + rhs_params
        #return '%s @= %s' % (lhs, rhs), params
        return '%s <@ %s AND %s @> %s' % (lhs, rhs, lhs, rhs), params + params

MolField.register_lookup(SameStructure)

################
# descriptors utils

def make_descriptor_mixin(name, prefix, field):
    return type(
        str('{0}_Mixin'.format(name.upper())),
        (object,),
        {
            'descriptor_name': name,
            'function': '{0}_{1}'.format(prefix, name),
            'default_output_field': field,
        },
    )


class DescriptorTransform(Transform):

    def as_sql(self, qn, connection):
        lhs, params = qn.compile(self.lhs)
        return "%s(%s)" % (self.function, lhs), params


##########################################
# MolField transforms and descriptors

MOL_DESCRIPTORS = [
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
    ('chi4v', FloatField),
    ('chi0n', FloatField),
    ('chi1n', FloatField),
    ('chi2n', FloatField),
    ('chi3n', FloatField),
    ('chi4n', FloatField),
    ('kappa1', FloatField),
    ('kappa2', FloatField),
    ('kappa3', FloatField),
    ('murckoscaffold', MolField),
]


MOL_DESCRIPTOR_MIXINS = [
    make_descriptor_mixin(d, 'mol', fieldkls())
    for d, fieldkls in MOL_DESCRIPTORS
]


MOL_DESCRIPTOR_TRANFORMS = [
    type(str('{0}_Transform'.format(mixin.descriptor_name.upper())),
         (mixin, DescriptorTransform,),
         { 'lookup_name': mixin.descriptor_name,
           'output_field': mixin.default_output_field, }
     )
    for mixin in MOL_DESCRIPTOR_MIXINS
]


for Transform in MOL_DESCRIPTOR_TRANFORMS:
    MolField.register_lookup(Transform)


##########################################
# RxnField transforms and descriptors

RXN_DESCRIPTORS = [
    ('numreactants', IntegerField),
    ('numproducts', IntegerField),
    ('numagents', IntegerField),
]


RXN_DESCRIPTOR_MIXINS = [
    make_descriptor_mixin(d, 'reaction', fieldkls())
    for d, fieldkls in RXN_DESCRIPTORS
]


RXN_DESCRIPTOR_TRANFORMS = [
    type(str('{0}_Transform'.format(mixin.descriptor_name.upper())),
         (mixin, DescriptorTransform,),
         { 'lookup_name': mixin.descriptor_name,
           'output_field': mixin.default_output_field, }
     )
    for mixin in RXN_DESCRIPTOR_MIXINS
]


for Transform in RXN_DESCRIPTOR_TRANFORMS:
    RxnField.register_lookup(Transform)


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
