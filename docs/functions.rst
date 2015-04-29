Database functions
==================

Most of the database functions implemented by the cartridge are also exposed through the django api and may be used to convert and create values, or in annotation and aggregate expressions::

  result = MoleculeModel.objects.aggregate(avg_amw=Avg(AMW('molecule')))

For consistency, all function names are defined as uppercase (this is not probably the prettiest solution, but it's easy to remember and unlikely to produce name clashes).

Functions
---------

- ``HBA``
- ``HBD``
- ``NUMATOMS``
- ``NUMHEAVYATOMS``
- ``NUMROTATABLEBONDS``
- ``NUMHETEROATOMS``
- ``NUMRINGS``
- ``NUMAROMATICRINGS``
- ``NUMALIPHATICRINGS``
- ``NUMSATURATEDRINGS``
- ``NUMAROMATICHETEROCYCLES``
- ``NUMAROMATICCARBOCYCLES``
- ``NUMALIPHATICCARBOCYCLES``
- ``NUMSATURATEDCARBOCYCLES``
- ``AMW``
- ``LOGP``
- ``TPSA``
- ``FRACTIONCSP3``
- ``CHI0V``
- ``CHI1V``
- ``CHI2V``
- ``CHI3V``
- ``CHI4V``
- ``CHI0N``
- ``CHI1N``
- ``CHI2N``
- ``CHI3N``
- ``CHI4N``
- ``KAPPA1``
- ``KAPPA2``
- ``KAPPA3``
- ``MURCKOSCAFFOLD``


- ``MOL``
- ``MOL_FROM_SMILES``
- ``MOL_FROM_SMARTS``
- ``MOL_FROM_CTAB``
- ``QMOL``
- ``QMOL_FROM_SMILES``
- ``QMOL_FROM_SMARTS``
- ``MOL_TO_SMILES``
- ``MOL_TO_SMARTS``
- ``MOL_TO_CTAB``
- ``MOL_INCHI``
- ``MOL_INCHIKEY``
- ``MOL_FORMULA``


- ``IS_VALID_SMILES``
- ``IS_VALID_SMARTS``
- ``IS_VALID_CTAB``


- ``NUMREACTANTS``
- ``NUMPRODUCTS``
- ``NUMAGENTS``


- ``REACTION``
- ``REACTION_FROM_SMILES``
- ``REACTION_FROM_SMARTS``
- ``REACTION_FROM_CTAB``
- ``REACTION_TO_SMILES``
- ``REACTION_TO_SMARTS``
- ``REACTION_TO_CTAB``
- ``REACTION_DIFFERENCE_FP``
- ``REACTION_STRUCTURAL_BFP``


- ``MORGAN_FP``
- ``MORGANBV_FP``
- ``FEATMORGAN_FP``
- ``FEATMORGANBV_FP``
- ``RDKIT_FP``
- ``ATOMPAIR_FP``
- ``ATOMPAIRBV_FP``
- ``TORSION_FP``
- ``TORSIONBV_FP``
- ``LAYERED_FP``
- ``MACCS_FP``


- ``TANIMOTO_SML``
- ``DICE_SML``
- ``TVERSKY_SML``
- ``TANIMOTO_DIST``
- ``DICE_DIST``

Aggregates
----------

- ``FMCS``
