SIESTA Bands workflow
++++++++++++++++++++++

Description
-----------

The **SiestaBandsWorkchain** workflow can be used to visualize the
electronic band structure of the system of interest. Its inputs are a
structure object and a specification of the quality and cost level of
the calculation. The latter is implemented internally, as in Quantum
Espresso, as a set of *protocols*, which group operational parameters
to offer the desired balance of accuracy and efficiency. Optionally,
the workflow will relax the geometry of the system before computing
the band structure. As discussed in the context of the Base workflow,
the computation could be implemented as a single (restartable) SIESTA
calculation, but it is instead segmented into different steps
(optional relaxation followed by a final electronic-structure plus
band calculation) to provide future support for different levels of
accuracy in the two stages. Support for the *fat-bands* feature that
tags energy levels with orbital projections will be added soon.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1 series, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).

Inputs
------

* **code**,  class :py:class:`Code  <aiida.orm.Code>`

A database object representing a Siesta executable.

* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`

A structure. See the plugin documentation for more details.

* **protocol**, class :py:class:`Str <aiida.orm.Str>`

Either "standard" or "fast" at this point.
Each has its own set of associated parameters.

- standard::

             {
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.2,
                'dm_convergence_threshold': 1.0e-4,
                'forces_convergence_threshold': "0.02 eV/Ang",
                'min_meshcutoff': 100, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'md-type-of-run': "cg",
                'md-num-cg-steps': 10,
                'pseudo_familyname': 'lda-ag',
                'atomic_heuristics': {
                    'H': { 'cutoff': 100 },
                    'Si': { 'cutoff': 100 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'DZP'
                }
	      }

- fast::
    
             {
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.25,
                'dm_convergence_threshold': 1.0e-3,
                'forces_convergence_threshold': "0.2 eV/Ang",
                'min_meshcutoff': 80, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'md-type-of-run': "cg",
                'md-num-cg-steps': 8,
                'pseudo_familyname': 'lda-ag',
                'atomic_heuristics': {
                    'H': { 'cutoff': 50 },
                    'Si': { 'cutoff': 50 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'SZP'
                }
	      }

The *atomic_heuristics* dictionary is intended to encode the
peculiarities of particular elements. It is work in progress.

The *basis* section applies globally for now.

Outputs
-------

* **output_parameters** :py:class:`Dict <aiida.orm.Dict>` 

A dictionary with metadata and scalar result values from the final *scf+bands*
calculation executed.

* **bands**, :py:class:`BandsData <aiida.orm.BandsData>`
  
Contains an array with the list of electronic energies for every
kpoint. For spin-polarized calculations, there is an extra dimension
for spin.

* **output_structure** :py:class:`StructureData <aiida.orm.StructureData>`
  
Present only if the workchain is modifying the geometry of the system.

