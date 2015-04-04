Model field reference
=====================

.. module:: django_rdkit.models.fields

.. currentmodule:: django_rdkit.models

A few custom model fields are provided, mapping to the chemical data types defined by the RDKit cartridge.

When fetched form the database, the python type of the model data attributes will match the corresponding RDKit class types, but for some fields additional data types may be used in assignment. For example, the value attributed to a ``MolField`` will be a ``Mol`` instance, and the value attributed to a ``ChemicalReaction`` will be a ``ChemicalReaction`` instance, but the values for these fields may be also assigned using a SMILES representation.

In trasferring data to and from the database, a binary representation is used for the fields supporting it.

Field types
-----------

``MolField``
............

.. class:: MolField(**options)

A field representing an RDKit molecule. It may be assigned using a ``Mol`` instance or with a SMILES string. Molecules values can be also created using one of the database functions implemented by the RDKit cartridge.


``RxnField``
............

.. class:: RxnField(**options)

A field storing a ``ChemicalReaction`` instance. It is assigned and returned from the database as a SMILES reactions string.


``BfpField``
............

.. class:: BfpField(**options)

A bit vector fingerprint. It may be assigned using an ``ExplicitBitVect`` instance or with an update query using one of the implemented fingerprint functions. 

``SfpField``
............

.. class:: SfpField(**options)


A sparse count vector fingerprint (``SparseIntVect``). Direct assignment of the ``SfpField`` field on the client side is not supported, the most practical way to assign values to the mapped table column is using an update query.







 
