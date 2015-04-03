Database setup
==============

The RDKit extension for PostgreSQL must be compiled and installed. 

Conda packages for the main RDKit releases, including the postgres cartridge for the linux platform, may be found from the `RDKit binstar channel <https://conda.binstar.org/rdkit>`_ or built using the recipes available from the `conda-rdkit repository <https://github.com/rdkit/conda-rdkit>`_.

Please refer to the `RDKit documentation <http://rdkit.readthedocs.org/en/latest/>`_ for general instructions regarding building the database cartridge from a source code distribution.

The details of the PostgreSQL database creation and configuration may vary depending on the deployment strategy and the application-specific needs, but no additional requirements exist for using the RDKit PostgreSQL cartridge in a django project. For general advice please refer to the official PostgreSQL and django documentation.

