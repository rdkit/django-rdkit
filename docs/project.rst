Project configuration and management
====================================

Compared to a regular django project using a postgresql database, using the RDKit cartridge mainly requires performing two additional operations.

Firstly, it is necessary to have the cartidge installed in the configured database. This operation corresponds to executing the following SQL statement::

  CREATE EXTENSION IF NOT EXISTS rdkit;

Moreover, efficient execution of structure and fingerprint based search operations requires the creation of an additional GiST index::

  CREATE INDEX index_name ON table_name USING GIST (column_name);

Two different mechanism are provided that allow integrating these operations with the development of the django application (in contrast to performing them independently, as part of the database administration).


The django-rdkit backend
------------------------

The ``django_rdkit`` package includes a custom database backend, implemented as a thin wrapper around the regular django ``postgresql_psycopg2`` backend. Usage is very simple and it just consists in configuring the database ``ENGINE`` settings as::

    DATABASES={
        "default": {
            'ENGINE': 'django_rdkit.backend',
	    ...
	}
    }


In preparing the database the custom backend will transparently ensure that the RDKit cartridge is configured. In addition, the database schema editor will recognize the model fields implemented by the ``django_rdkit`` package and it will automatically include the creation of a corresponding GiST index everytime a chemical field is used (unless the field ``chem_index`` option is set as ``False``).


Migration Operations
--------------------

.. module:: django_rdkit.operations

.. currentmodule:: django_rdkit.operations

Alternatively, it is possible to configure the project using the usual ``postgresql_psycopg2`` backend, and perform the same administration tasks with the help of some custom operations in the database migrations. In particular, two operations are implemented:


.. class:: RDKitExtension( )

An ``Operation`` subclass that will install the RDKit cartridge.


.. class:: GiSTIndex(model_name, name, index_name=None)

An ``Operation`` sublcass that wraps the management of a GiST index for field ``name`` of model ``model_name`` (please note that when the custom ``django_rdkit`` backend is not used the value of the ``chem_index`` field option has no effect). The optional ``index_name`` parameter allows customizing the name of the created index.

Compared to using the ``django_rdkit`` custom backend described above, migration operations might provide some additional flexibility and a more precise control on when changes are applied to the schema (for example, the backend will add a GiST index only within the context of a model creation or field addition). 
