from __future__ import unicode_literals

from django.db import models

INTEGER_DESCRIPTORS = [
    'hba', 
    'hbd',
    'numatoms',
    'numheavyatoms',
    'numrotatablebonds',
    'numheteroatoms',
    'numrings',
    'numaromaticrings',
    'numaliphaticrings',
    'numsaturatedrings',
    'numaromaticheterocycles',
    'numaliphaticheterocycles',
    'numsaturatedheterocycles',
    'numaromaticcarbocycles',
    'numaliphaticcarbocycles',
    'numsaturatedcarbocycles',
]


FLOAT_DESCRIPTORS = [
    'amw',
    'logp',
    'tpsa',
    'fractioncsp3',
    'chi0v', 'chi1v', 'chi2v', 'chi3v', 'chi5v',
    'chi0n', 'chi1n', 'chi2n', 'chi3n', 'chi5n',
    'kappa1', 'kappa2', 'kappa3', 'kappa4',
]


def _make_mixin(name, field):
    return type(
        str('{0}_Mixin'.format(name.upper())), 
        (object,),
        { 
            'descriptor_name': name, 
            'function': 'mol_{0}'.format(name),
            'output_field': field, 
        },
    )


DESCRIPTOR_MIXINS = (
    [ _make_mixin(d, models.IntegerField()) for d in INTEGER_DESCRIPTORS ] +
    [ _make_mixin(d, models.FloatField()) for d in FLOAT_DESCRIPTORS ]
    )
 
