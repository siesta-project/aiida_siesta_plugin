# Changelog

## v1.4.0

Version compatible with aiida-core>=1.3.0,<2.0.0.

### Improvements
- Introduce several workchains to help the optimization of basis sets.

### For developers
- Introduces dependences on aiida-optimize.
- It is supposed to be the last release before moving the development to 2.0.x (that will be compatible
  with aiida-core>=2.0.0 only).

## v1.3.0

Version compatible with aiida-core>=1.3.0,<2.0.0.

### Improvements
- Introduce dedicated support for optical calculations. In particular,
  dedicated inputs are created for the `SiestaCalculation` and `SiestaBaseWorkChain`.
  Moreover a `EpsilonWorkChain` workchain is created to automatically
  obtain the electronic contribution to the static dielectric constant.
- Introduce `SiestaBaseNEBWorkChain`, the core building block for the creation of workflows that enable
  the search of the Minimum Energy Pathway (MEP) connecting two local minima of the potential energy surface
  through the Nudge Elastic Band (NEB) method.

### Bug fixes
- Clarify units of stress.
- The `SiestaBaseWorkChain` now correcty handle LUA inputs.

## v1.2.0

Version compatible with aiida-core>=1.3.0,<2.0.0.
Support for python 3.6 has been dropped.

### Improvements
- Enable LUA interface to SIESTA. New inputs have been introduce to activate LUA control
  over siesta calculations.
- Access to the [sisl](https://github.com/zerothi/sisl) functionalities, since sisl is now
  dependency of this package.
- Introduction of the `IonData` data type that host the entity that in siesta terminology is called “ion”.
  An "ion" packages the set of basis orbitals and KB projectors for a given species and therefore
  can be use in substitution of the pseudo and basis inputs.
  An `IonData` node is also returned in output for each species of the calculation so to help
  the improve the tracking of used orbitals and KB projectors.
  From an `IonData` instance, the list of orbitals of the calculation can be easily accessed in the form
  of `SislAtomicOrbitals`.
- Migration to the use of the [aiida-pseudo](https://github.com/aiidateam/aiida-pseudo) package for the
  management of pseudo data types and families. The development on the the support for pseudopotentials
  types has been centralized in an external package. The use of aiida-siesta custom classes `PsmlData`
  and `PsfData` has been deprecated and will be remove in version 2.0.0.

### Bug fixes
- Fixed the atom heuristic selection of basis in the protocol system

### For developers
- Heavy refactoring to take advantage of the aiida-core validators system.
- Start reorganization of folders. Especially moved the location of the protocols
  system and of tkdict features.
- Updated dependencies of tests and pre-commits.

## v1.1.1

Version compatible with aiida-core>=1.3.0,<2.0.0.

### Improvements
- Improvement in the `BandGapWorkChain`. A new feature has been introduced:
  when no `bandskpoints` are set in input, the workchain automatically calculates the bands
  using a k-space path automatically selected by SeeK-path.
  In particular, if a molecular dynamics run was also requested (usually a relaxation),
  the dynamics is run first and the bands are calculated on a separate Siesta run. This
  allows to select the k-space path for bands automatically using the output structure.
  When the `bandskpoints` is selected in input, the behavior is the same of previous versions.
- Added an output of the `SiestaSequentialConverger` listing, if any,
  the unconverged parameters. More explicit warnings are returned when a parameter
  does not converge.
- A tutorial is now available in the documentation, covering the submission of simple siesta
  calculations, the use of protocols, the use of the iterator/converger workchains and
  a simple example on how to create custom workchains.

### Bug fixes
- A bug was reported regarding the wrong selection of the kpoints path for bands
  when a relaxation with variable cell was also requested. The energy values were correct,
  but the distance between kpoints was not corresponding to the one of the Siesta
  file ".bands". The bug is now fixed.
- The method `inputs_generator().get_inputs_dict` (available for `SiestaCalculations` and
  all the workchains in the package) had a minor bug. In case the kpoints or bandskpoints
  were not requested, the dictionary was anyway returning items for "bandskpoints"
  and "kpoints" with value "None". This has been fixed.
  The method `inputs_generator().get_filled_builder` was instead correct.
- A bug was leading to the failure of `SiestaSequentialConverger` in case one of
  the parameters to converge in the sequential process was not reaching convergence.
  It is now fixed and, in addition, an output listing the unconverged parameters has been
  added.

### For developers
- Set up a nightly built in GitHub actions that tests the plugin (run all its tests)
  against the develop branch of `aiida-core`.
- The `BandGapWorkChain` has been refactored almost entirely.
- The abstract class `SequentialConverger` has been modified to fix the bug
  of `SiestaSequentialConverger` described above.

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
