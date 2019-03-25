SIESTA STM workflow
++++++++++++++++++++++

Description
-----------

The **SiestaVibraWorkchain** workflow produces files with the phonon
frequencies and eigenvectors from a previous Siesta calculation.

The inputs to the vibra workchain include the Siesta code, the Vibra
code, a (possibly relaxed) structure, a supercell structure, a protocol
and a reciprocal space path for plotting the phonon bands.

The Vibra package included in the Siesta release contains a series of
programs to compute the phonon frequencies: 1) FCBuild, which builds
a supercell to compute the dynamical matrix and 2) Vibra, which reads
the force constants and computes the phonon spectrum. The first part
(FCBuild) is substituted in the workflow by a direct calculation of
the supercell in the input file (see ``examples/workflows/test_vibra_wk.py``).
After the supercell is created and stored as a 
:py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>` object,
the workflow starts. First it launches a Siesta calculation to calculate
the force constants and then launches a Vibra calculation to compute
the phonon spectrum.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1 series, which
can be found in the development platform
(http://launchpad.net/siesta/).

Inputs
------

* **code**, a code associated to the Siesta plugin

* **vibra_code**, a code associated to the Vibra plugin

* **structure**, class :py:class:`StructureData
  <aiida.orm.data.structure.StructureData>`

A structure. See the plugin documentation for more details.

* **structure_sc**, class :py:class:`StructureData
  <aiida.orm.data.structure.StructureData>`

A super-cell structure built from the previous structure.

* **protocol**, Str

Either "standard" or "fast" at this point.
Each has its own set of associated parameters.

- standard::

             {
                'kpoints_mesh': 1,
                'dm_convergence_threshold': 1.0e-4,
                'min_meshcutoff': 100, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'pseudo_familyname': 'si_ldapsf',
                'atomic_heuristics': {
                    'Si': { 'cutoff': 100 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'DZP'
                }
	      }

- fast::
    
             {
                'kpoints_mesh': 1,
                'dm_convergence_threshold': 1.0e-3,
                'min_meshcutoff': 80, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'pseudo_familyname': 'si_ldapsf',
                'atomic_heuristics': {
                    'Si': { 'cutoff': 50 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'SZ'
                }
	      }

The *basis* section applies globally for now.

Outputs
-------

Two files, one with the phonon frequencies (aiida.bands) and the other
with the phonon eigenvectors (aiida.vectors).



