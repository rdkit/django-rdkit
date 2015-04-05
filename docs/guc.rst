Configuration parameters
========================

The RDKit PostgreSQL cartridge defines a set of configuration parameters fine-tuning some of the implemented functions. From the SQL interface these parameters can be manipulated using ``set`` and ``show`` statements::

  my_database=# set rdkit.do_chiral_sss=true;

but in the ``django_rdkit`` package they are also exposed as attributes of a ``config`` object::

  In [43]: from django_rdkit.config import config

so that their values can be set and queried without leaving the python domain::

  In [44]: config.do_chiral_sss = True

  In [45]: print(config.tanimoto_threshold)
  0.5

As you may notice from the examples, the main difference should consist in a minor change in the naming convention. The cartridge defines these parameters with a name starting with an ``rdkit.`` prefix. In naming the corresponding attributes of the ``config`` object this prefix is dropped.
