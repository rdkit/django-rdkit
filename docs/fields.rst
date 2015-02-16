Model field reference
=====================

.. module:: django_rdkit.db.models.fields

.. currentmodule:: django_rdkit.db.models

A few custom model fields are provided, mapping to the chemical data types defined by the RDKit cartridge.

Field options
-------------

In addition to those usually supported by django fields, the following field options are supported:

``chem_index``
..............

If ``True``, a GIST index will be created for this column. Default value is ``True``.


Field types
-----------

``MolField``
............

.. class:: MolField(**options)

A field representing an RDKit molecule. It may be assigned using a ``Mol`` instance or with a SMILES string. Molecules values can be also created using one of the database functions implemented by the RDKit cartridge.


``BfpField``
............

.. class:: BfpField(**options)

A bit vector fingerprint. It may be assigned using an ``ExplicitBitVect`` instance or with an update query using one of the implemented fingerprint functions. 

``SfpField``
............

.. class:: SfpField(**options)


A sparse count vector fingerprint (``SparseIntVect``). Values for this field are assigned with an update query using one of the implemented fingerprint functions.



 
