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

    def test_insert_null(self):
        MoleculeModel.objects.create(molecule=None)
        objs = MoleculeModel.objects.filter(molecule=None)
        self.assertEqual(objs.count(), 1)

    def test_assign_from_sql_expression(self):
        mol = MoleculeModel.objects.create(molecule=MOL_FROM_SMILES(Value('c1ccccc1')))
        updated = MoleculeModel.objects.filter(id=mol.id).update(molecule=MOL_FROM_SMILES(Value('c1ccccc1')))
        self.assertEqual(updated, 1)

    def test_exact_lookup(self):

        objs = MoleculeModel.objects.filter(molecule='COC(c1ccccc1)c1ccccc1')
        self.assertEqual(objs.count(), 1)

        objs = MoleculeModel.objects.filter(molecule='Nc1ccc(Cl)nc1')
        self.assertEqual(objs.count(), 1)

        objs = MoleculeModel.objects.filter(molecule=Chem.MolFromSmiles('Nc1ccc(Cl)nc1'))
        self.assertEqual(objs.count(), 1)

        objs = MoleculeModel.objects.filter(molecule=MOL_FROM_SMILES(Value('Nc1ccc(Cl)nc1')))
        self.assertEqual(objs.count(), 1)

    def test_hassubstruct_lookup(self):

        objs = MoleculeModel.objects.filter(molecule__hassubstruct='C1=C(C)C=CC=C1')
        self.assertEqual(objs.count(), 61)

        objs = MoleculeModel.objects.filter(
            molecule__hassubstruct=MOL_FROM_SMILES(Value('C1=C(C)C=CC=C1')))
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

        objs = MoleculeModel.objects.filter(molecule__hassubstruct=Chem.MolFromSmiles('C1=CN=CC=C1'))
        cnt4 = objs.count()
        self.assertEqual(cnt2, 7)

    def test_issubstruct_lookup(self):

        objs = MoleculeModel.objects.filter(molecule__issubstruct='CCN1c2ccccc2Sc2ccccc21')
        self.assertEqual(objs.count(), 2)

        objs = MoleculeModel.objects.filter(molecule__issubstruct='CC[N+]([O-])(CC)CCCN1c2ccccc2S(=O)c2ccccc21')
        self.assertEqual(objs.count(), 4)

        objs = MoleculeModel.objects.filter(molecule__issubstruct=Chem.MolFromSmiles('CC[N+]([O-])(CC)CCCN1c2ccccc2S(=O)c2ccccc21'))
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

        aggr = MoleculeModel.objects.aggregate(sum_amw=Sum(AMW('molecule')))
        self.assertAlmostEqual(aggr['sum_amw'], 23687.4, 3)

        aggr = MoleculeModel.objects.aggregate(max_amw=Max(AMW('molecule')))
        self.assertAlmostEqual(aggr['max_amw'], 836.468, 3)

        aggr = MoleculeModel.objects.aggregate(min_amw=Min(AMW('molecule')))
        self.assertAlmostEqual(aggr['min_amw'], 94.497, 3)

    def test_exercise_mol_descriptor_transforms(self):
        # this is a test that doesn't test much, it simply executes queries,
        # hoping for regressions
        for desc in (
                'hba',
                'hbd',
                'numatoms',
                'numheavyatoms',
                'numrotatablebonds',
                'numheteroatoms',
                'numrings',
                'numaromaticrings',
                'numaliphaticrings',
                'numsaturatedrings',
                'numaromaticheterocycles',
                'numaliphaticheterocycles',
                'numsaturatedheterocycles',
                'numaromaticcarbocycles',
                'numaliphaticcarbocycles',
                'numsaturatedcarbocycles',
        ):
            kwargs = {
                'molecule__{}__gt'.format(desc): 0,
            }
            _ = MoleculeModel.objects.filter(**kwargs).count()

        for desc in (
                'amw',
                'logp',
                'tpsa',
                'fractioncsp3',
                'chi0v',
                'chi1v',
                'chi2v',
                'chi3v',
                'chi4v',
                'chi0n',
                'chi1n',
                'chi2n',
                'chi3n',
                'chi4n',
                'kappa1',
                'kappa2',
                'kappa3',
        ):
            kwargs = {
                'molecule__{}__gt'.format(desc): 0.,
            }
            _ = MoleculeModel.objects.filter(**kwargs).count()

        kwargs = {
            'molecule__murckoscaffold__amw__gt': 50.
        }
        _ = MoleculeModel.objects.filter(**kwargs).count()

    def test_exercise_mol_functions(self):
        for func in (
                HBA,
                HBD,
                NUMATOMS,
                NUMHEAVYATOMS,
                NUMROTATABLEBONDS,
                NUMHETEROATOMS,
                NUMRINGS,
                NUMAROMATICRINGS,
                NUMALIPHATICRINGS,
                NUMSATURATEDRINGS,
                NUMAROMATICHETEROCYCLES,
                NUMAROMATICCARBOCYCLES,
                NUMALIPHATICCARBOCYCLES,
                NUMSATURATEDCARBOCYCLES,

                AMW,
                LOGP,
                TPSA,
                FRACTIONCSP3,
                CHI0V,
                CHI1V,
                CHI2V,
                CHI3V,
                CHI4V,
                CHI0N,
                CHI1N,
                CHI2N,
                CHI3N,
                CHI4N,
                KAPPA1,
                KAPPA2,
                KAPPA3,
                MURCKOSCAFFOLD,

                MOL_TO_SMILES,
                MOL_TO_SMARTS,
                MOL_TO_CTAB,
                MOL_INCHI,
                MOL_INCHIKEY,
                MOL_FORMULA,
                MOL_TO_SVG,
        ):
            _ = list(MoleculeModel.objects.annotate(result=func('molecule')))


class RxnFieldTest(TestCase):

    def setUp(self):
        for smiles in REACTION_SMILES_SAMPLE:
            ReactionModel.objects.create(rxn=smiles)
        for smarts in REACTION_SMARTS_SAMPLE:
            ReactionModel.objects.create(
                rxn=Chem.ReactionFromSmarts(str(smarts)))

    def test_products_count(self):
        self.assertEqual(
            ReactionModel.objects.filter(rxn__numproducts=1).count(),
            6
        )
        self.assertEqual(
            ReactionModel.objects.filter(rxn__numproducts__gt=1).count(),
            2
        )

    def test_reactants_count(self):
        self.assertEqual(
            ReactionModel.objects.filter(rxn__numreactants=2).count(),
            2
        )
        self.assertEqual(
            ReactionModel.objects.filter(rxn__numreactants__lt=2).count(),
            6
        )

    def test_agents_count(self):
        self.assertEqual(
            ReactionModel.objects.filter(rxn__numagents=0).count(),
            4
        )
        self.assertEqual(
            ReactionModel.objects.filter(rxn__numagents__gt=1).count(),
            2
        )

    def test_descriptors(self):
        qs = ReactionModel.objects.all()

        aggr = qs.aggregate(sum_prods=Sum(NUMPRODUCTS('rxn')))
        self.assertEqual(aggr['sum_prods'], 10)

        aggr = qs.aggregate(max_prods=Max(NUMPRODUCTS('rxn')))
        self.assertEqual(aggr['max_prods'], 2)

        aggr = qs.aggregate(min_prods=Min(NUMPRODUCTS('rxn')))
        self.assertEqual(aggr['min_prods'], 1)

    def test_exercise_rxn_descriptor_transforms(self):
        # another test that doesn't test much, it simply executes queries,
        # hoping for regressions
        for desc in (
                'numreactants',
                'numproducts',
                'numagents',
        ):
            kwargs = {
                'rxn__{}__gt'.format(desc): 0,
            }
            _ = ReactionModel.objects.filter(**kwargs).count()

    def test_exercise_rxn_functions(self):
        for func in (
                NUMREACTANTS,
                NUMPRODUCTS,
                NUMAGENTS,
                REACTION_TO_SMILES,
                REACTION_TO_SMARTS,
                REACTION_TO_CTAB,
        ):
            _ = list(ReactionModel.objects.annotate(result=func('rxn')))


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


class IsValidSmilesTest(TestCase):

    def setUp(self):
        for smiles in ('c1ccccc1', 'whatever', 'c1cccc1', 'CCCO'):
            _ = SmilesModel.objects.create(smiles=smiles)

    def test_is_valid_smiles(self):
        qs = SmilesModel.objects.annotate(valid=IS_VALID_SMILES('smiles'))
        valid = qs.values_list('valid', flat=True)
        self.assertTrue(valid[0])
        self.assertFalse(valid[1])
        self.assertFalse(valid[2])
        self.assertTrue(valid[3])


class IsValidCtabTest(TestCase):

    def setUp(self):
        mol = Chem.MolFromSmiles('c1cocc1')
        CtabModel.objects.create(ctab=Chem.MolToMolBlock(mol))
        CtabModel.objects.create(ctab='rubbish')

    def test_is_valid_ctab(self):
        qs = CtabModel.objects.annotate(valid=IS_VALID_CTAB('ctab'))
        valid = qs.values_list('valid', flat=True)
        self.assertTrue(valid[0])
        self.assertFalse(valid[1])


class IsValidSmartsTest(TestCase):

    def setUp(self):
        for smarts in ('c1ccc[c,n]c1',
                       'whatever',
                       '(F)F.[c1:1][c:2rubbish',
                       'C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]'):
            _ = SmartsModel.objects.create(smarts=smarts)

    def test_is_valid_smiles(self):
        qs = SmartsModel.objects.annotate(valid=IS_VALID_SMARTS('smarts'))
        valid = qs.values_list('valid', flat=True)
        self.assertTrue(valid[0])
        self.assertFalse(valid[1])
        self.assertFalse(valid[2])
        self.assertTrue(valid[3])


class MolFromSmilesTest(TestCase):

    def setUp(self):
        for smiles in ('c1ccccc1', 'whatever', 'c1cccc1', 'CCCO',):
            _ = SmilesModel.objects.create(smiles=smiles)

    def test_mol_from_smiles(self):
        qs = SmilesModel.objects.annotate(valid=IS_VALID_SMILES('smiles'))
        qs = qs.filter(valid=True)
        updated = qs.update(molecule=MOL_FROM_SMILES('smiles'))
        self.assertEqual(updated, 2)
        isnull_count = SmilesModel.objects.filter(molecule__isnull=True).count()
        self.assertEqual(isnull_count, 2)
