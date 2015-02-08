
$ createdb django_rdkit_tutorial

$ django-admin startproject tutorial

$ cd tutorial/

$ emacs tutorial/settings.py &

DATABASES={
    "default": {
        'ENGINE': 'django_rdkit.db.backends.postgresql_psycopg2',
        'NAME': 'django_rdkit_tutorial',
        'USER': '',
        'PASSWORD': '',
        'HOST': 'localhost',
        'PORT': '',
    }
}

$ python manage.py migrate
$ python manage.py runserver

$ python manage.py startapp app

$ emacs app/models.py

$ emacs tutorial/settings.py &

$ python manage.py makemigrations app

$ python manage.py sqlmigrate app 0001

$ python manage.py check

$ python manage.py migrate

$ python manage.py shell

In [1]: path = '../../chembl/chembl_20_chemreps.txt'

In [2]: from rdkit import Chem

In [3]: def chembl(path):
   ...:     with open(path, 'rt') as f:
   ...:         f.next() # skip header
   ...:         for line in f:
   ...:             tokens = line.split()
   ...:             name, smiles = tokens[:2]
   ...:             molecule = Chem.MolFromSmiles(smiles)
   ...:             if molecule:
   ...:                 yield name, molecule
   ...:                 

In [4]: from app.models import Compound

In [7]: Compound.objects.bulk_create(
(Compound(molecule=molecule) for molecule in molecules(path)),
batch_size=10000
)

In [6]: for molecule in molecules(path):
   ...:     Compound.objects.create(molecule=molecule)
   ...:    

In [16]: for name, molecule in chembl(path):
    smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
    m = Chem.MolFromSmiles(smiles)
    if not m:
        print "interconversion issue for", name
        continue
    Compound.objects.create(molecule=smiles)

In [19]: print Compound.objects.count()
1455711

In [24]: Compound.objects.filter(molecule__hassubstruct='c1cccc2c1nncc2').count()
Out[24]: 392

In [25]: Compound.objects.filter(molecule__hassubstruct='c1ccnc2c1nccn2').count()
Out[25]: 782

In [26]: Compound.objects.filter(molecule__hassubstruct='c1cncc2n1ccn2').count() 
Out[26]: 1234

In [27]: Compound.objects.filter(molecule__hassubstruct='Nc1ncnc(N)n1').count()
Out[27]: 5395

In [28]: Compound.objects.filter(molecule__hassubstruct='c1scnn1').count()
Out[28]: 12665

In [29]: Compound.objects.filter(molecule__hassubstruct='c1cccc2c1ncs2').count() 
Out[29]: 15977

In [30]: Compound.objects.filter(molecule__hassubstruct='c1cccc2c1CNCCN2').count()
Out[30]: 1512

In [33]: for c in Compound.objects.filter(molecule__hassubstruct='c1cccc2c1CNCCN2')[:100]:
    print c.name, Chem.MolToSmiles(c.molecule)
   ....:     
CHEMBL3276818 O=C1c2ccccc2N(c2ccccc2)CCN1CCN1CCCCC1
CHEMBL3277031 Cc1ccc2c(c1)C(=O)NCC(=O)N2C
CHEMBL3277027 CN1c2c(cc(Cl)cc2[N+](=O)[O-])C(=O)NCC1=O
CHEMBL3277028 CN1c2ccc(Cl)cc2C(=O)NCC1=O
CHEMBL3277029 CN1c2cc(Cl)ccc2C(=O)NCC1=O
...
CHEMBL3134101 CCC1(C)NC(=O)c2cc(S(=O)(=O)Nc3ccc(C(F)(F)F)cc3OC)cc(Cl)c2NC1=O
CHEMBL3134102 CCC1(C)NC(=O)c2cc(S(=O)(=O)Nc3ccc(C)cc3C#N)cc(Cl)c2NC1=O
CHEMBL3134103 CC1(C)NC(=O)c2c(F)c(S(=O)(=O)Nc3ccc(C(F)(F)F)cc3)ccc2NC1=O
CHEMBL3134104 CC1(C)NC(=O)c2cc(S(=O)(=O)Nc3ccc(C(F)(F)F)cc3)cc(F)c2NC1=O
CHEMBL3134105 CC1(C)NC(=O)c2cc(S(=O)(=O)Nc3ccc(C(F)(F)F)cc3)c(F)cc2NC1=O

In [34]: from django_rdkit.db.models import QMOL, Value

In [39]: query = QMOL(Value('c1[o,s]ncn1'))

In [38]: for c in Compound.objects.filter(molecule__hassubstruct=query)[:500]:    print c.name, Chem.MolToSmiles(c.molecule)
   ....:     
CHEMBL602115 Cc1nc(-c2c(Cl)cc(Cl)cc2-c2cnc(C(C)NC(=O)N(C)O)c(F)c2)no1
CHEMBL1082636 O=C(c1ccco1)N1CSCC1c1nc(-c2cccc(Cl)c2)no1
CHEMBL1076357 COc1ccccc1CNC(=O)c1cc(C(F)(F)F)nn1-c1cccc(-c2noc(C(C)N)n2)c1
CHEMBL563301 COC(=O)CCc1nc(C2CC(c3ccc(O)c(F)c3)=NO2)no1
CHEMBL7083 CC(=CCn1oc(=O)[nH]c1=O)c1cccc(OCc2noc(-c3ccccc3)n2)c1
...
CHEMBL105541 Cc1cccc(-c2noc(CN3C(=O)c4ccccc4C3=O)n2)c1
CHEMBL108856 Cc1nc(-c2ccc(CSc3nc4ccccc4n3Cc3ccc(Cl)cc3)cc2)no1
CHEMBL107705 Cc1noc(C2CN(C)CCC2c2ccc(Cl)cc2)n1
CHEMBL108598 Cc1noc(C2CNCCC2c2ccc(Cl)cc2)n1
CHEMBL115982 O=c1[nH]c(=O)n(C(CCc2ccc(OCc3cccc(Cl)c3)cc2)c2ccccc2)o1

# No use of stereochemistry yet

# Similarity searches will use some fingerprints

$ emacs app/models.py

$ python manage.py makemigrations

$ python manage.py migrate

$ python manage.py shell

In [4]: from django_rdkit.db.models import *

In [6]: from app.models import Compound

In [3]: Compound.objects.update(
   ...: torsionbv=TORSIONBV_FP('molecule'),
   ...: mfp2=MORGANBV_FP('molecule'),
   ...: ffp2=FEATMORGANBV_FP('molecule'),
   ...: )


