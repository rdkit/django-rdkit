from django.db import connection

_PARAMETERS = (
    'tanimoto_threshold',
    'dice_threshold',
    'do_chiral_sss',
    'do_enhanced_stereo_sss',
    'sss_fp_size',
    'morgan_fp_size',
    'featmorgan_fp_size',
    'layered_fp_size',
    'rdkit_fp_size',
    'hashed_torsion_fp_size',
    'hashed_atompair_fp_size',
)

class Config(object):

    def __getattr__(self, name):
        if not name in _PARAMETERS:
            raise AttributeError
        with connection.cursor() as c:
            c.execute("SHOW rdkit.{}".format(name))
            value = c.fetchone()[0]
        return value

    def __setattr__(self, name, value):
        if not name in _PARAMETERS:
            raise AttributeError
        with connection.cursor() as c:
            c.execute("SET rdkit.{}=%s".format(name), (value,))


config = Config()
