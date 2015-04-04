Project configuration and management
====================================

Django projects integrating the functionalities provided by ``django_rdkit`` should configure their settings to use the PostgreSQL database backend::

  DATABASES={
      'default': {
          'ENGINE': 'django.db.backends.postgresql_psycopg2',
          # [...],
      }
  }

Two additional operations are then to be performed on the database in order to use the RDKit cartridge.

Firstly, it is necessary to have the cartidge installed in the configured database. This operation corresponds to executing the following SQL statement::

  CREATE EXTENSION IF NOT EXISTS rdkit;

One simple way to integrate this operation within a django project consists in installing the ``django_rdkit`` package as a django application::

  INSTALLED_APPS = (
      # [...]
      'django_rdkit',
  )

When the database is configured a migration will be this way executed, ensuring that that the RDKit extension is created.

Moreover, efficient execution of structure and similarity searches requires the creation of an additional GiST index::

  CREATE INDEX index_name ON table_name USING GIST (column_name);

The creation of this custom index is also supported with a migration operation. Users should include this operation in the implementation of a custom migration for the fields that may require it.


Migration Operations
--------------------
.. module:: django_rdkit.operations

.. currentmodule:: django_rdkit.operations

.. class:: RDKitExtension( )

An ``Operation`` subclass that will install the RDKit cartridge (install ``django_rdkit`` as a django application to include a migration that wraps this operation).


.. class:: GiSTIndex(model_name, name, index_name=None)

An ``Operation`` subclass that wraps the management of a GiST index for field ``name`` of model ``model_name``. The optional ``index_name`` parameter allows customizing the name of the created index.
