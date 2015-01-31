#from __future__ import unicode_literals

from django.core.management import call_command
from django.test import TestCase

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

        objs = TestModel.objects.filter(molecule__hassubstruct='C1=CN=CC=C1')
        self.assertEqual(objs.count(), 7)

    def test_issubstruct_lookup(self):

        objs = TestModel.objects.filter(molecule__issubstruct='CCN1c2ccccc2Sc2ccccc21')
        self.assertEqual(objs.count(), 2)

        objs = TestModel.objects.filter(molecule__issubstruct='CC[N+]([O-])(CC)CCCN1c2ccccc2S(=O)c2ccccc21')
        self.assertEqual(objs.count(), 4)

        
        
