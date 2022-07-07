## [0.3.2]

### Fixed
- The `MolField` and related lookup expressions were refactored to support
  bulk update operations

## [0.3.1]

### Fixed
- Replaced the deprecated/unsupported `ugettext_lazy` with `gettext_lazy`

## [0.3.0]

### Added
- `mol_to_svg` support

### Removed
- The gist index migration operation was removed (GiST support is available in django)

## [0.2.0]

### Changed
- Updated for Django 3.0
- Support for Python2 was removed

## [Unreleased]

### Added
- README and CHANGELOG files

### Changed
- MolField values are always passed to the database as pickled Mol instances
- MolField supports parsing mol blocks and InChI input in addition to SMILES
- MolField specifies an appropriate default form field input

## [0.1.0] - 2018-08-16
