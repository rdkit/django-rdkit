from __future__ import unicode_literals

from django.test import TestCase

from django_rdkit.models import *

from rdkit.Chem import AllChem as Chem

from .models import *
from .molecules import SMILES_SAMPLE
from .reactions import REACTION_SMILES_SAMPLE, REACTION_SMARTS_SAMPLE

class MolFieldTest(TestCase):
    
    def setUp(self):
        for smiles in SMILES_SAMPLE:
            MoleculeModel.objects.create(molecule=smiles)

    def test_exact_lookup(self):

        objs = MoleculeModel.objects.filter(molecule='COC(c1ccccc1)c1ccccc1')
        self.assertEqual(objs.count(), 1)

        objs = MoleculeModel.objects.filter(molecule='Nc1ccc(Cl)nc1')
        self.assertEqual(objs.count(), 1)

    def test_hassubstruct_lookup(self):

        objs = MoleculeModel.objects.filter(molecule__hassubstruct='C1=C(C)C=CC=C1')
        self.assertEqual(objs.count(), 61)

        objs = MoleculeModel.objects.filter(molecule__hassubstruct='C1=CC=CC=C1')
        cnt1 = objs.count()
        self.assertEqual(cnt1, 70)

        objs = MoleculeModel.objects.filter(molecule__hassubstruct='C1=CN=CC=C1')
        cnt2 = objs.count()
        self.assertEqual(cnt2, 7)

        objs = MoleculeModel.objects.filter(
            Q(molecule__hassubstruct='C1=CC=CC=C1') |
            Q(molecule__hassubstruct='C1=CN=CC=C1'),
        )
        cnt3 = objs.count()
        self.assertEqual(cnt3, 73)
        self.assertTrue(cnt3 <= cnt1 + cnt2)

        qmol = QMOL(Value('c1[c,n]cccc1'))
        objs = MoleculeModel.objects.filter(molecule__hassubstruct=qmol)
        self.assertEqual(objs.count(), cnt3)

    def test_issubstruct_lookup(self):

        objs = MoleculeModel.objects.filter(molecule__issubstruct='CCN1c2ccccc2Sc2ccccc21')
        self.assertEqual(objs.count(), 2)

        objs = MoleculeModel.objects.filter(molecule__issubstruct='CC[N+]([O-])(CC)CCCN1c2ccccc2S(=O)c2ccccc21')
        self.assertEqual(objs.count(), 4)

    def test_descriptor_AMW(self):

        threshold = 250.
        cnt1 = MoleculeModel.objects.filter(molecule__amw__gt=threshold).count()
        self.assertEqual(cnt1, 36)
        
        objs = MoleculeModel.objects.annotate(amw=AMW('molecule'))
        cnt2 = sum(1 for m in objs if m.amw > threshold)
        self.assertEqual(cnt2 - cnt1, 0)

        aggr = MoleculeModel.objects.aggregate(avg_amw=Avg(AMW('molecule')))
        self.assertAlmostEqual(aggr['avg_amw'], 236.874, 3)


class RxnFieldTest(TestCase):

    def setUp(self):
        for smiles in REACTION_SMILES_SAMPLE:
            ReactionModel.objects.create(rxn=smiles)
        for smarts in REACTION_SMARTS_SAMPLE:
            ReactionModel.objects.create(rxn=Chem.ReactionFromSmarts(smarts))

    def test_placeholder(self):
        self.assertEqual(
            ReactionModel.objects.count(), 
            len(REACTION_SMILES_SAMPLE) + len(REACTION_SMARTS_SAMPLE)
        )


class BfpFieldTest1(TestCase):
    
    def setUp(self):
        for smiles in SMILES_SAMPLE:
            record = BfpModel.objects.create()
            record.bfp = MORGANBV_FP(Value(smiles))
            record.save()

    def test_tanimoto_lookup(self):
        query_bfp = MORGANBV_FP(Value('CCN1c2ccccc2Sc2ccccc21'))
        objs = BfpModel.objects.filter(bfp__tanimoto=query_bfp)
        self.assertEqual(objs.count(), 2)

    def test_dice_lookup(self):
        query_bfp = MORGANBV_FP(Value('CCN1c2ccccc2Sc2ccccc21'))
        objs = BfpModel.objects.filter(bfp__dice=query_bfp)
        self.assertEqual(objs.count(), 5)

        
class BfpFieldTest2(TestCase):
    
    def setUp(self):
        for smiles in SMILES_SAMPLE:
            record = BfpModel.objects.create()
            record.bfp=MORGANBV_FP(Value(smiles), Value(3))
            record.save()

    def test_tanimoto_lookup(self):
        query_bfp = MORGANBV_FP(Value('CCN1c2ccccc2Sc2ccccc21'), Value(3))
        objs = BfpModel.objects.filter(bfp__tanimoto=query_bfp)
        self.assertEqual(objs.count(), 2)

    def test_dice_lookup(self):
        query_bfp = MORGANBV_FP(Value('CCN1c2ccccc2Sc2ccccc21'), Value(3))
        objs = BfpModel.objects.filter(bfp__dice=query_bfp)
        self.assertEqual(objs.count(), 2)


class BfpFieldTest3(TestCase):

    def test_pkl_io(self):
        bfps = {}
        for smiles in SMILES_SAMPLE:
            mol = Chem.MolFromSmiles(smiles)
            bfp = Chem.GetMorganFingerprintAsBitVect(mol, 2, 512)
            obj = BfpModel.objects.create(bfp=bfp)
            bfps[obj.pk] = bfp

        for obj in BfpModel.objects.all():
            self.assertTrue(obj.pk in bfps)
            ibfp = bfps[obj.pk]
            obfp = obj.bfp
            self.assertEqual(list(ibfp.GetOnBits()),
                             list(obfp.GetOnBits()))


class SfpFieldTest1(TestCase):
    
    def setUp(self):
        for smiles in SMILES_SAMPLE:
            record = SfpModel.objects.create()
            record.sfp = MORGAN_FP(Value(smiles))
            record.save()

    def test_tanimoto_lookup(self):
        query_sfp = MORGAN_FP(Value('CCN1c2ccccc2Sc2ccccc21'))
        objs = SfpModel.objects.filter(sfp__tanimoto=query_sfp)
        self.assertEqual(objs.count(), 2)

    def test_dice_lookup(self):
        query_sfp = MORGAN_FP(Value('CCN1c2ccccc2Sc2ccccc21'))
        objs = SfpModel.objects.filter(sfp__dice=query_sfp)
        self.assertEqual(objs.count(), 15)

        
class SfpFieldTest2(TestCase):
    
    def setUp(self):
        for smiles in SMILES_SAMPLE:
            record = SfpModel.objects.create()
            record.sfp = MORGAN_FP(Value(smiles), Value(5))
            record.save()

    def test_tanimoto_lookup(self):
        query_sfp = MORGAN_FP(Value('CCN1c2ccccc2Sc2ccccc21'), Value(5))
        objs = SfpModel.objects.filter(sfp__tanimoto=query_sfp)
        self.assertEqual(objs.count(), 2)

    def test_dice_lookup(self):
        query_sfp = MORGAN_FP(Value('CCN1c2ccccc2Sc2ccccc21'), Value(5))
        objs = SfpModel.objects.filter(sfp__dice=query_sfp)
        self.assertEqual(objs.count(), 2)


        
