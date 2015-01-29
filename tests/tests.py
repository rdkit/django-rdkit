from django.core.management import call_command
from django.test import TestCase

from .models import TestModel

class MigrationsTests(TestCase):
    def test_makemigrations(self):
        call_command('makemigrations', dry_run=True)


class MoleculeFieldTest(TestCase):
    
    def setUp(self):

        molecules = [
            'c1ccccc1',
            'c1ccccc1C',
            'c1cc(C)ccc1C',
            'n1ccccc1',
            'n1ccccc1C',
            'n1cc(C)ccc1C',
            'n1ccccc1',
            'C1CCCCC1',
        ]

        for molecule in molecules:
            TestModel.objects.create(molecule=molecule)

    def test_exact_lookup(self):

        count = TestModel.objects.filter(molecule='C1=CC=CC=C1').count()
        self.assertEqual(count, 1)

        count = TestModel.objects.filter(molecule='C1C=CCCC1').count()
        self.assertEqual(count, 0)

        count = TestModel.objects.filter(molecule='C1=CN=CC=C1').count()
        self.assertEqual(count, 2)

    def test_hassubstruct_lookup(self):

        objs = TestModel.objects.filter(molecule__hassubstruct='C1=C(C)C=CC=C1')
        self.assertEqual(objs.count(), 2)

        objs = TestModel.objects.filter(molecule__hassubstruct='C1=CN=CC=C1')
        self.assertEqual(objs.count(), 4)

    def test_issubstruct_lookup(self):

        objs = TestModel.objects.filter(molecule__issubstruct='C1=C(C)C=CC=C1')
        self.assertEqual(objs.count(), 2)

        objs = TestModel.objects.filter(molecule__issubstruct='C1=C(C)N=CC=C1')
        self.assertEqual(objs.count(), 3)


        
        
