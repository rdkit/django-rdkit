Field Lookups
=============

MolField
--------

Lookup operators
................

- ``hassubstruct``
- ``issubstruct``
- ``exact``

Descriptor transforms
.....................

Most of the molecular descriptor functions defined by the cartridge are also
available as transform operators. To ease the mnemonics, the name of these operators is based on the original function name, deprived of the ``mol_`` prefix (``mol_hba`` becomes ``hba``) and following the usual django conventions all names are lowercase. For example::

  # count all compounds with AMW above a provided threshold value
  CompoundModel.objects.filter(molecule__amw__gt=threshold).count()


- ``hba``
- ``hbd``
- ``numatoms``
- ``numheavyatoms``
- ``numrotatablebonds``
- ``numheteroatoms``
- ``numrings``
- ``numaromaticrings``
- ``numaliphaticrings``
- ``numsaturatedrings``
- ``numaromaticheterocycles``
- ``numaliphaticheterocycles``
- ``numsaturatedheterocycles``
- ``numaromaticcarbocycles``
- ``numaliphaticcarbocycles``
- ``numsaturatedcarbocycles``
- ``amw``
- ``logp``
- ``tpsa``
- ``fractioncsp3``
- ``chi0v``
- ``chi1v``
- ``chi2v``
- ``chi3v``
- ``chi4v``
- ``chi0n``
- ``chi1n``
- ``chi2n``
- ``chi3n``
- ``chi4n``
- ``kappa1``
- ``kappa2``
- ``kappa3``
- ``murckoscaffold``


RxnField
--------

Lookup operators
................

- ``hassubstruct``
- ``hassubstructfp``
- ``issubstruct``
- ``issubstructfp``

Descriptor transforms
.....................

- ``numreactants``
- ``numproducts``
- ``numagents``

BfpField and SfpField
---------------------

Lookup operators
................

- ``tanimoto``
- ``dice``


 
