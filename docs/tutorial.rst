Tutorial
========

This tutorial will try to reproduce the operations described in the `RDKit PostgreSQL cartridge documentation <http://rdkit.readthedocs.org/en/latest/Cartridge.html>`_, but within the context of a django project.

Some familiarity with django and the django database api is assumed (excellent documentation about these is available from the `django web site <https://docs.djangoproject.com>`_).

PostgreSQL and the RDKit cartridge should be installed and running on the system. A database should be created with appropriate access privileges to be used by the tutorial project. Minimally, this requires running the following command::

  $ createdb django_rdkit_tutorial


Creation of the tutorial project
--------------------------------

Create a new skeleton django project named ``tutorial_project``::

  $ django-admin startproject tutorial_project

Change working directory to the ``tutorial_project`` directory (where the ``manage.py`` file is located) and open the ``tutorial_project/settings.py`` module with your favourite text editor. 

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

And extend the ``INSTALLED_APPS`` list to include the ``django_rdkit`` application::

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
  
The ``migrate`` command above configures the database for the installed applications. The inclusion of ``django_rdkit`` in the ``INSTALLED_APP`` is not strictly required, but allows integrating the creation of the RDKit extension with the management of the django project, as evidenced by using ``sqlmigrate``::

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

The additional functionalities developed in the context of this tutorial will be contained in a so-called django application. We'll call this application ``tutorial_application``::

  $ python manage.py startapp tutorial_application

The list of ``INSTALLED_APPS`` in the ``tutorial_project/settings.py`` module must be extended to include the new application::

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

We'll use this application to manage a collection of compound structures. In order to do so, edit the ``tutorial_applications/models.py`` module so that it looks like the following::

  from django_rdkit import models
  
  class Compound(models.Model):
  
      name = models.CharField(max_length=256)
      molecule = models.MolField()

Please note that we import ``models`` from the ``django_rdkit`` package, instead of from ``django.db`` as we would usually do. This makes the ``MolField`` and the other functionalities that are specific the RDKit cartridge available, together with the rest of the usual fields and functions that are usually availble from ``django.db``.

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

  
Structures import and substructure queries
------------------------------------------

To display the use of structure searches we'll use a copy of the ChEMBL data. Download a copy of the ``chembl_20_chemreps.txt`` which is available from `here <ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_20/>`_ and place it into a suitable directory.

The initial import may therefore be performed with code similar to the following::

  $ python manage.py shell
  [...]
  In [1]: path = '../../chembl/chembl_20_chemreps.txt'
   
  In [2]: from rdkit import Chem
  
  In [3]: def chembl(path, limit=None):
     ...:     count = 0
     ...:     with open(path, 'rt') as f:
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

The import loop may take some time, consider using the ``limit`` parameter to shorten the duration of this step. Once the import has completed one can easily verify the number of available compounds::

  In [8]: Compound.objects.count()
  Out[8]: 1455712

In order to efficiently perform structural queries on the imported compounds, a database index must be created. This operation can be implemented with adding a special internal class ``Meta`` to the model and assigning the ``indexes`` the corresponding index::

  from django_rdkit import models
  from django.contrib.postgres.indexes import GistIndex

  class Compound(models.Model):
  
      name = models.CharField(max_length=256)
      molecule = models.MolField()
      
      class Meta:
        indexes = [
            GistIndex(fields=['molecule']),
        ]

Execute the following command to create a migration::

  $ python manage.py makemigrations --name create_compound_molecule_index tutorial_application
  Migrations for 'tutorial_application':
    0002_create_compound_molecule_index.py:

When done, save your changes and run the migration (depending on the number of structures imported into the model, the indexing may take quite some time to complete)::

  $ python manage.py migrate tutorial_application
  Operations to perform:
    Apply all migrations: tutorial_application
  Running migrations:
    Rendering model states... DONE
    Applying tutorial_application.0002_create_compound_molecule_index...


Finally, following the original tutorial, we can now perform a few example substructure queries::

  In [1]: from django_rdkit.models import *
  
  In [2]: from tutorial_application.models import *
  
  In [3]: def smiles_substructure_query(substructure):
     ....:     query = Compound.objects.filter(molecule__hassubstruct=substructure)
     ....:     for cmpd in query.annotate(smiles=MOL_TO_SMILES('molecule'))[:5]:
     ....:         print(cmpd.name, cmpd.smiles)
     ....:         

The above code uses the ``hassubstruct`` lookup operator, which is specific to the ``MolField`` field, and also uses the ``MOL_TO_SMILES`` database function to convert the selected molecules and annotate the model instance with a smiles string. Both functionalities are provided by the RDKit cartridge.

::

  In [4]: smiles_substructure_query('c1cccc2c1nncc2')
  CHEMBL113970 CCCCn1c(=O)c2cc(OC)c(OC)cc2c2nnc3cc4c(cc3c21)OCO4
  CHEMBL113470 COc1cc2c(cc1OC)c1nnc3cc4c(cc3c1n(C(C)CN(C)C)c2=O)OCO4
  CHEMBL12112 CC(C)Sc1ccc(CC2CCN(C3CCN(C(=O)c4cnnc5ccccc54)CC3)CC2)cc1
  CHEMBL71086 COc1cc2c(cc1OC)c1nnc3cc4c(cc3c1n(CCN(C)C)c2=O)OCO4 
  CHEMBL89981 c1ccc(CN2CCC(CCNc3cc4ccc5ccccc5c4nn3)CC2)cc1
  
  In [5]: smiles_substructure_query('c1ccnc2c1nccn2')
  CHEMBL110168 CCOC(=O)Nc1cc(NC(C)CCCN(CC)CC)c2nc(-c3ccccc3)c(-c3ccccc3)nc2n1
  CHEMBL50456 Clc1ccc(CN2CCN(c3nc4cccnc4n4cccc34)CC2)c(Cl)c1
  CHEMBL107535 O=c1c2cccn2c2ncccc2n1CCNC(=S)Nc1ccc(Br)cn1
  CHEMBL51225 c1cc2c(N3CCN(c4ccccc4)CC3)nc3cccnc3n2c1
  CHEMBL54246 Cc1ccnc2c1nc(N1CCN(Cc3ccccc3)CC1)c1cccn12

SMARTS-based queries
....................

Similarly, substructure queries can use a SMARTS string as argument::

  In [20]: def smarts_substructure_query(substructure):
     ....:     query = Compound.objects.filter(molecule__hassubstruct=QMOL(Value(substructure)))
     ....:     for cmpd in query.annotate(smiles=MOL_TO_SMILES('molecule'))[:5]:
     ....:         print(cmpd.name, cmpd.smiles)
     ....:         

The lookup api expects a SMILES string by default, so a query molecule must be created explicitly, using the ``QMOL`` constructor, which is exposed as a database function. Please note that database functions execute on the backend, and by default assume their argument to resolve to a database column. Since a literal SMARTS string is used, it must be wrapped inside a call to ``Value()`` (the query expression api was introduced in django 1.8, for further details about this see the official `documentation <https://docs.djangoproject.com/en/1.8/ref/models/expressions/>`_.

::

  In [21]: smarts_substructure_query('c1[o,s]ncn1')
  CHEMBL52013 C[C@@H](NC(=O)c1nsc(-c2ccc(Cl)cc2)n1)[C@](O)(Cn1cncn1)c1ccc(F)cc1F
  CHEMBL48759 CCN(CC)C(=O)N1Cc2c(-c3noc(C4CC4)n3)ncn2-c2ccccc21
  CHEMBL48839 CCSC(=O)N1Cc2c(-c3noc(C4CC4)n3)ncn2-c2ccccc21
  CHEMBL105111 COc1ccc(-c2noc(CN3C(=O)c4ccccc4C3=O)n2)cc1
  CHEMBL105112 Cc1ccccc1-c1noc(CN2C(=O)c3ccccc3C2=O)n1

Using stereochemistry
.....................

By default stereochemistry is not taken into account when performing substructure queries::

  In [42]: smiles_substructure_query('NC(=O)[C@H]1CCCN1C=O')
  CHEMBL118176 CC(C)[C@@H](NC(=O)COc1ccc(OCC(=O)O)cc1)C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)c1nc2ccccc2o1)C(C)C
  CHEMBL117981 O=C(CCCc1ccccc1)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)c1ccccn1
  CHEMBL117920 O=C(CCCc1ccccc1)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)c1cccnc1
  CHEMBL117024 Cc1ccc(C[C@H](NC(=O)c2ccc(C)c(O)c2C)[C@H](O)C(=O)N2C[C@@H](Cl)C[C@H]2C(=O)NC(C)(C)C)cc1
  CHEMBL117088 Cc1cc(O)c(C)c(C(=O)N[C@@H](Cc2cccc(C(F)(F)F)c2)[C@H](O)C(=O)N2C[C@@H](Cl)C[C@H]2C(=O)NC(C)(C)C)c1

As described in the RDKit documentation, the cartridge defines a set of configuration parameters that allow controlling this and other aspects. These parameters are exposed as attributes of a ``config`` object::

  In [43]: from django_rdkit.config import config


In particular, the effect of stereochemistry on the results returned by substructure searches is changed using the ``do_chiral_sss`` configuration variable::

  In [44]: config.do_chiral_sss = True
  
  In [45]: smiles_substructure_query('NC(=O)[C@H]1CCCN1C=O')
  CHEMBL100712 N=C(N)NCCC[C@H]1NC(=O)[C@H]2CCCN2C(=O)[C@H](Cc2ccccc2)NC(=O)CCCCCCCCCCNC(=O)C1=O
  CHEMBL98474 Cc1ccccc1S(=O)(=O)NC(=O)N1CCC[C@@H]1C(=O)NCCC(=O)NC(Cc1c[nH]cn1)C(=O)O
  CHEMBL2369135 CC[C@H](C)[C@@H]1NC(=O)[C@@H]([C@H](C)c2c(C)cc(OC)cc2C)NC(=O)[C@H](N)C(C)(C)SSC[C@@H]2NC(=O)[C@@H](CC(N)=O)NC(=O)[C@@H](CCC(=O)NCCCC[C@H](C(=O)NCC(N)=O)NC(=O)[C@H]3CCCN3C2=O)NC1=O
  CHEMBL2369136 CC[C@H](C)[C@@H]1NC(=O)[C@H]([C@@H](C)c2c(C)cc(OC)cc2C)NC(=O)[C@H](N)C(C)(C)SSC[C@@H]2NC(=O)[C@@H](CC(N)=O)NC(=O)[C@@H](CCC(=O)NCCCC[C@H](C(=O)NCC(N)=O)NC(=O)[C@H]3CCCN3C2=O)NC1=O
  CHEMBL98856 N=C(N)NCCC[C@H]1NC(=O)[C@H]2CCCN2C(=O)[C@H](Cc2ccccc2)NC(=O)CCCCCNC(=O)C1=O


Similarity queries
------------------

Open the file ``tutorial_application/models.py`` for editing again, and extend the ``Compound`` model with some fingerprint fields, as displayed below::

  from django_rdkit import models
  
  class Compound(models.Model):
  
      name = models.CharField(max_length=256)
      molecule = models.MolField()
  
      torsionbv = models.BfpField(null=True)
      mfp2 = models.BfpField(null=True)
      ffp2 = models.BfpField(null=True)

(please note that the new fields are defined as nullable so that we can alter the existing database table adding initially empty columns).

Create a corresponding schema migration::

  $ python manage.py makemigrations tutorial_application --name add_compound_fingerprint_fields
  Migrations for 'tutorial_application':
    0003_add_compound_fingerprint_fields.py:
      - Add field ffp2 to compound
      - Add field mfp2 to compound
      - Add field torsionbv to compound

And finally, apply it to the current schema::

  $ python manage.py migrate tutorial_application
  Operations to perform:
    Apply all migrations: tutorial_application
  Running migrations:
    Rendering model states... DONE
    Applying tutorial_application.0003_add_compound_fingerprint_fields... 

The fingerpring columns may be filled with data that is computed with an update query::

  $ python manage.py shell
  [...]
  In [1]: from django_rdkit.models import *
  
  In [2]: from tutorial_application.models import Compound 
  
  In [3]: Compound.objects.update(
     ...: torsionbv=TORSIONBV_FP('molecule'),
     ...: mfp2=MORGANBV_FP('molecule'),
     ...: ffp2=FEATMORGANBV_FP('molecule'),
     ...: )
  Out[3]: 1455712

Once this query has completed, an index must still be added on the column (or columns) that will be frequently used to perform similarity queries. This database administration step may be again be another entry in the ``indexes`` field of the ``Meta`` internal class. First, add an index to the ``indexes``::

  from django_rdkit import models
  from django.contrib.postgres.indexes import GistIndex
  
  class Compound(models.Model):
  
      name = models.CharField(max_length=256)
      molecule = models.MolField()
      
      torsionbv = models.BfpField(null=True)
      mfp2 = models.BfpField(null=True)
      ffp2 = models.BfpField(null=True)

      class Meta:
        indexes = [
            GistIndex(fields=['molecule']),
            GistIndex(fields=['mfp2']),
        ]

Second, create a migration::

  $ python manage.py makemigrations --name create_compound_mfp2_index tutorial_application
  Migrations for 'tutorial_application':
    0004_create_compound_mfp2_index.py:

And then run the migration to complete the preparation of the database::

  $ python manage.py migrate tutorial_application
  Operations to perform:
    Apply all migrations: tutorial_application
  Running migrations:
    Rendering model states... DONE
    Applying tutorial_application.0004_create_compound_mfp2_index...

The following demonstrate a basic similarity search::

  In [1]: from django_rdkit.models import *
  
  In [2]: from tutorial_application.models import *
  
  In [3]: smiles = 'Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1'
  
  In [4]: value = MORGANBV_FP(Value(smiles))
  
  In [5]: Compound.objects.filter(mfp2__tanimoto=value).count()Out[6]: 67

Following the original tutorial from the RDKit documentation, the next step consists in implementing a query to return the sorted list of neighbors along with the accompanying SMILES::

  In [8]: def get_mfp2_neighbors(smiles):
     ...:     value = MORGANBV_FP(Value(smiles))
     ...:     queryset = Compound.objects.filter(mfp2__tanimoto=value)
     ...:     queryset = queryset.annotate(smiles=MOL_TO_SMILES('molecule'))
     ...:     queryset = queryset.annotate(sml=TANIMOTO_SML('mfp2', value))
     ...:     queryset = queryset.order_by(TANIMOTO_DIST('mfp2', value)) 
     ...:     queryset = queryset.values_list('name', 'smiles', 'sml')
     ...:     return queryset
     ...: 

The function wraps a non-trivial database api expression, but the generated SQL query can be easily displayed for a sample queryset::

  In [22]: qs = get_mfp2_neighbors('c1ccccc1')
  
  In [23]: print(qs.query)
  SELECT "tutorial_application_compound"."name",
  mol_to_smiles("tutorial_application_compound"."molecule") AS "smiles",
  tanimoto_sml("tutorial_application_compound"."mfp2", morganbv_fp(c1ccccc1)) AS "sml"
  FROM "tutorial_application_compound" WHERE
  "tutorial_application_compound"."mfp2" % (morganbv_fp(c1ccccc1)) ORDER BY
  ("tutorial_application_compound"."mfp2" <%> morganbv_fp(c1ccccc1)) ASC

You can use the ``get_mfp2_neighbors`` function to perform some sample queries::

  In [9]: for name, smiles, sml in get_mfp2_neighbors('Cc1ccc2nc(-c3ccc(NC(C4N(C(c5cccs5)=O)CCC4)=O)cc3)sc2c1')[:10]:
      print(name, smiles, sml)
     ...:     
  CHEMBL467428 Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(C(=O)c5cccs5)CC4)cc3)sc2c1 0.772727272727273
  CHEMBL461435 Cc1ccc2nc(-c3ccc(NC(=O)C4CCCN(S(=O)(=O)c5cccs5)C4)cc3)sc2c1 0.657534246575342
  CHEMBL460340 Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(S(=O)(=O)c5cccs5)CC4)cc3)sc2c1 0.647887323943662
  CHEMBL460588 Cc1ccc2nc(-c3ccc(NC(=O)C4CCN(S(=O)(=O)c5cccs5)C4)cc3)sc2c1 0.638888888888889
  CHEMBL1608585 O=C(Nc1nc2ccc(Cl)cc2s1)[C@@H]1CCCN1C(=O)c1cccs1 0.623188405797101
  CHEMBL1327784 COc1ccc2nc(NC(=O)[C@@H]3CCCN3C(=O)c3cccs3)sc2c1 0.619718309859155
  CHEMBL518028 Cc1ccc2nc(-c3ccc(NC(=O)C4CN(S(=O)(=O)c5cccs5)C4)cc3)sc2c1 0.611111111111111
  CHEMBL1316870 Cc1ccc(NC(=O)C2CCCN2C(=O)c2cccs2)cc1C 0.606060606060606
  CHEMBL1309021 O=C(Nc1ccc(S(=O)(=O)N2CCCC2)cc1)C1CCCN1C(=O)c1cccs1 0.602941176470588
  CHEMBL1706764 Cc1ccc(NC(=O)C2CCCN2C(=O)c2cccs2)c(C)c1 0.597014925373134
  
  In [10]: for name, smiles, sml in get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1')[:10]:
     ....:     print(name, smiles, sml)
     ....:     
  CHEMBL394654 Cc1ccc2nc(N(C)CCN(C)c3nc4ccc(C)cc4s3)sc2c1 0.692307692307692
  CHEMBL491074 CN(CC(=O)O)c1nc2cc([N+](=O)[O-])ccc2s1 0.583333333333333
  CHEMBL1617304 CC(=O)N(CCCN(C)C)c1nc2ccc(C)cc2s1 0.571428571428571
  CHEMBL1350062 CC(=O)N(CCCN(C)C)c1nc2ccc(C)cc2s1.Cl 0.549019607843137
  CHEMBL1621941 Cc1ccc2nc(N(CCN(C)C)C(=O)c3cc(Cl)sc3Cl)sc2c1 0.518518518518518
  CHEMBL1626442 Cc1ccc2nc(N(CCCN(C)C)C(=O)CS(=O)(=O)c3ccccc3)sc2c1 0.517857142857143
  CHEMBL1617545 Cc1ccc2nc(N(CCCN(C)C)C(=O)CCc3ccccc3)sc2c1 0.517857142857143
  CHEMBL406760 Cc1ccc2nc(NC(=O)CCC(=O)O)sc2c1 0.510204081632653
  CHEMBL1624740 Cc1ccc(S(=O)(=O)CC(=O)N(CCCN(C)C)c2nc3ccc(C)cc3s2)cc1 0.509090909090909
  CHEMBL1620007 Cc1ccc2nc(N(CCN(C)C)C(=O)c3ccc4ccccc4c3)sc2c1 0.509090909090909


Adjusting the similarity cutoff
...............................

::

  In [11]: print(get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1').count())
  18
  
  In [12]: from django_rdkit.config import config
  
  In [13]: config.tanimoto_threshold = 0.7
  
  In [14]: print(get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1').count())
  0
  
  In [15]: config.tanimoto_threshold = 0.6
  
  In [16]: print(get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1').count())
  1
  
  In [17]: config.tanimoto_threshold = 0.5
  
  In [18]: print(get_mfp2_neighbors('Cc1ccc2nc(N(C)CC(=O)O)sc2c1').count())
  18


