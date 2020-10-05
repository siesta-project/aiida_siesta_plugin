# Changelog

## v1.1.0

Version compatible with aiida-core>=1.3.0,<2.0.0. Previous versions of aiida-core do not have a working `BaseRestartWorkChain`, which is the building piece of many workchain in aiida-siesta.
Support for python 2.7 and python 3.5 has been dropped.

### Improvements
- Extend the `STMCalculation` to work with modern versions of the "plstm" code. Included support for spin options.
- Improved the Siesta parser to recognize more siesta errors.
- Refactoring of the `SiestaBaseWorkChain` to make use of `BaseRestartWorkChain` from aiida-core, added a more complete errors handling.
- Substitute the bands workchain with the `BandGapWorkChain`. No parameter is assumed inside the workchian ("protocols" are moved outside).
- Change the `STMWorkChain`. No parameter is assumed inside the workchain and it has been extended to support all the modes and spins options for STM calculations.
- Introduce several new workchains: `SiestaIterator`, `SiestaConverger`, `SiestaSequentialConverger` and `EqOfStateFixedCellShape`.
- Introduce the "protocols" system, that suggests inputs for AiiDA workchains based on the structure under investigation and a "protocol" string.
- Add several examples about new features, workchains and the protocol system.
- Expanded the documentation.

### For developers
- Introduce a complete set of unitests.
- Introduced proper CI with GitHub actions, including pre-commit requirements.
- Improved the `FDFDict`: more methods and a behavior more similar to a python dictionary.
- Refactoring of `PsmlData` and `PsfData` classes to use new Group feature of aiida-core.


## v1.0.1

Version compatible with aiida-core>=1.0.0,<2.0.0

### Improvements
- `More robust SiestaBaseWorkChain`.

### For developers
- Avoid the re-setting of Siesta filepaths via spec.inp.
- Fix obsolete call and input attributes in `SiestaBaseWorkChain`.
- Remove obsolete files in several places.
- Raise exception in parser if XML file is faulty.

## v1.0.0

This version is the first release of aiida-siesta that is compatible with aiida-core 1.0.0. Previous versions of aiida-core are not compatible with this release.
