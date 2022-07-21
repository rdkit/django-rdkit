Django settings
===============

The following django settings are supported:

``DJANGO_RDKIT_MOL_SERIALIZATION``
..................................

Default: ``'BINARY'``

Controls whether the value associated to a ``MolField`` is transferred to/from the database as a binary serialized (pickled) RDKit molecule, or as a SMILES text string. Supported values are: ``BINARY``, ``TEXT``.
