Database setup
==============

The RDKit extension for PostgreSQL must be installed. Please refer to the `RDKit wiki <http://code.google.com/p/rdkit/wiki/BuildingTheCartridge>`_ for detailed instructions.

creation of a chemical database template
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**NOTE: This is no longer necessary. The django-rdkit backend will automatically install the extension during the database preparation step.**

::

    # Creating the template database.
    $ createdb -E UTF8 template_rdkit
    # Allows non-superusers the ability to create from this template
    $ psql -d postgres -c "UPDATE pg_database SET datistemplate='true' WHERE datname='template_rdkit';"
    # Loading the RDKit cartridge
    $ psql -d template_rdkit -c "CREATE EXTENSION rdkit"

creation of the database
^^^^^^^^^^^^^^^^^^^^^^^^

::

    $ createdb -T template_rdkit my_chem_db

