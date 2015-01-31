#from __future__ import unicode_literals

from django.core.management import call_command
from django.test import TestCase

from django_rdkit.db.models import Q, Value, QMOL

from .models import TestModel
from .molecules import SMILES_SAMPLE

class MoleculeFieldTest(TestCase):
    
    def setUp(self):
        for smiles in SMILES_SAMPLE:
            TestModel.objects.create(molecule=smiles)

    def test_exact_lookup(self):

        objs = TestModel.objects.filter(molecule='COC(c1ccccc1)c1ccccc1')
        self.assertEqual(objs.count(), 1)

        objs = TestModel.objects.filter(molecule='Nc1ccc(Cl)nc1')
        self.assertEqual(objs.count(), 1)

    def test_hassubstruct_lookup(self):

        objs = TestModel.objects.filter(molecule__hassubstruct='C1=C(C)C=CC=C1')
        self.assertEqual(objs.count(), 61)

        objs = TestModel.objects.filter(molecule__hassubstruct='C1=CC=CC=C1')
        cnt1 = objs.count()
        self.assertEqual(cnt1, 70)

        objs = TestModel.objects.filter(molecule__hassubstruct='C1=CN=CC=C1')
        cnt2 = objs.count()
        self.assertEqual(cnt2, 7)

        objs = TestModel.objects.filter(
            Q(molecule__hassubstruct='C1=CC=CC=C1') |
            Q(molecule__hassubstruct='C1=CN=CC=C1'),
        )
        cnt3 = objs.count()
        self.assertEqual(cnt3, 73)
        self.assertTrue(cnt3 <= cnt1 + cnt2)

        qmol = QMOL(Value('c1[c,n]cccc1'))
        objs = TestModel.objects.filter(molecule__hassubstruct=qmol)
        self.assertEqual(objs.count(), cnt3)

    def test_issubstruct_lookup(self):

        objs = TestModel.objects.filter(molecule__issubstruct='CCN1c2ccccc2Sc2ccccc21')
        self.assertEqual(objs.count(), 2)

        objs = TestModel.objects.filter(molecule__issubstruct='CC[N+]([O-])(CC)CCCN1c2ccccc2S(=O)c2ccccc21')
        self.assertEqual(objs.count(), 4)

        
        
