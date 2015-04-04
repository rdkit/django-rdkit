Tutorial
========

This tutorial will try to reproduce the operations described in the `RDKit PostgreSQL cartridge documentation <http://rdkit.readthedocs.org/en/latest/Cartridge.html>`_, but within the context of a django project.

Some familiarity with django and the django database api is assumed (excellent documentation about these is available from the `django web site <https://docs.djangoproject.com>`_).

Database setup
--------------

PostgreSQL and the RDKit cartridge should be installed and running on the system. 

A database should be created with appropriate access privileges to be used by the tutorial project. Minimally, this requires running the following command::

  $ createdb django_rdkit_tutorial


Creation of the django project
------------------------------

Select an appropriate filesystem location and create a new skeleton django project named `tutorial_project`::

  $ django-admin startproject tutorial_project

Change working directory to the `tutorial_project` directory (where the `manage.py` file is located) and open the `tutorial_project/settings.py` module with your favourite text editor. 

Replace the default database settings with those appropriate to the created PostgreSQL database::

  DATABASES={
      'default': {
          'ENGINE': 'django.db.backends.postgresql_psycopg2',
          'NAME': 'django_rdkit_tutorial',
          'USER': '',
          'PASSWORD': '',
          'HOST': '',
          'PORT': '',
      }
  }

And extend the `INSTALLED_APPS` list to include the `django_rdkit` application::

  INSTALLED_APPS = (
      'django.contrib.admin',
      'django.contrib.auth',
      'django.contrib.contenttypes',
      'django.contrib.sessions',
      'django.contrib.messages',
      'django.contrib.staticfiles',
      'django_rdkit',
  )


Finally, initialize the database::

  $ python manage.py migrate
  
The `migrate` command above configures the database for the installed applications. The inclusion of `django_rdkit` in the `INSTALLED_APP` is not strictly required, but allows integrating the creation of the RDKit extension with the management of the django project, as evidenced by using `sqlmigrate`::

  $ python manage.py sqlmigrate django_rdkit 0001
  BEGIN;
  CREATE EXTENSION IF NOT EXISTS rdkit;
  
  COMMIT;
 
The correct configuration of the project may be quickly verified at this stage by running a direct SQL query using the database connection that is created by django::

  $ python manage.py shell
  [...]
  In [1]: from django.db import connection
  
  In [2]: with connection.cursor() as cursor:
     ...:     cursor.execute("SELECT mol_amw('C')")
     ...:     print(cursor.fetchone()[0])
     ...:     
  16.043


Creation of a django application
--------------------------------

The additional functionalities developed in the context of this tutorial will be contained in a so-called django application. We'll call this application `tutorial_application'::

  $ python manage.py startapp tutorial_application

The list of `INSTALLED_APPS` in the `tutorial_project/settings.py` module must be extended to include the new application::

  INSTALLED_APPS = (
      'django.contrib.admin',
      'django.contrib.auth',
      'django.contrib.contenttypes',
      'django.contrib.sessions',
      'django.contrib.messages',
      'django.contrib.staticfiles',
      'django_rdkit',
      'tutorial_application',
  )

We'll use this application to manage a collection of compound structures. In order to do so, edit the `tutorial_applications/models.py` module so that it looks like the following::

  from django_rdkit import models
  
  class Compound(models.Model):
  
      name = models.CharField(max_length=256)
      molecule = models.MolField()

Please note that we import `models` from the `django_rdkit` package, instead of from `django.db` as we would usually do. This makes the `MolField` and the other functionalities that are specific the RDKit cartridge available, together with the rest of the usual fields and functions that are usually availble from `django.db`.

In order to extend the schema of the PostgreSQL database to include this model, we now need to create and apply a corresponding migration::

  $ python manage.py makemigrations tutorial_application
  Migrations for 'tutorial_application':
    0001_initial.py:
      - Create model Compound
  $ python manage.py migrate tutorial_application
  Operations to perform:
    Apply all migrations: tutorial_application
  Running migrations:
    Rendering model states... DONE
    Applying tutorial_application.0001_initial... OK

We can immediately try adding data to this model using again the python shell::

  $ python manage.py shell
  [...]
  In [1]: from tutorial_application.models import Compound
  
  In [2]: Compound.objects.create(name='benzene', molecule='c1ccccc1')
  Out[2]: <Compound: Compound object>
  
  In [3]: from django_rdkit.models import *
  
  In [4]: for compound in Compound.objects.annotate(amw=AMW('molecule')):
     ...:     print(compound.name, compound.amw)
     ...:     
  benzene 78.114

We can now delete this sample compound, more data will be imported in the next section of this tutorial::

  In [5]: Compound.objects.all().delete()

  
Structures import and search
----------------------------

To display the use of structure searches we'll use a copy of the ChEMBL data. Download a copy of the `chembl_20_chemreps.txt` which is available from `here <ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_20/>`_ and place it into a suitable directory.

The initial import may therefore be performed with code similar to the following::

  $ python manage.py shell
  [...]
  In [1]: path = '../../chembl/chembl_20_chemreps.txt'
   
  In [2]: from rdkit import Chem
  
  In [3]: def chembl(path, limit=None):
     ...:     count = 0
     ...:     with open(path, 'rt) as f:
     ...:         next(f) # skip header
     ...:         for line in f:
     ...:             name, smiles = line.split()[:2]
     ...:             molecule = Chem.MolFromSmiles(smiles)
     ...:             if molecule:
     ...:                 yield name, molecule
     ...:                 count += 1
     ...:             if limit and count == limit:
     ...:                 break
     ...:             
  
  In [4]: from tutorial_application.models import Compound
  
  In [5]: for name, molecule in chembl(path, limit=None): 
     ...:     smiles = Chem.MolToSmiles(molecule)
     ...:     test_molecule = Chem.MolFromSmiles(smiles)
     ...:     if not test_molecule:
     ...:         print('smiles-mol-smiles roundtrip issue:', name)
     ...:     else:
     ...:         Compound.objects.create(name=name, molecule=molecule)
     ...:         

The import loop may take some time, consider using the `limit` parameter to shorten the duration of this step.

