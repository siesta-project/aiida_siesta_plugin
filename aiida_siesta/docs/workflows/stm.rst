SIESTA STM workflow
++++++++++++++++++++++

Description
-----------

The **SiestaSTMWorkchain** workflow is functionally very similar
to the **SiestaBandsWorkchain** workflow, but instead of a band
structure, the analysis stage produces a file with the local density
of states (LDOS) in an energy window. The LDOS can be seen as a
"partial charge density" to which only those wavefunctions with
eigenvalues in a given energy interval contribute. In the
Tersoff-Hamann approximation, the LDOS can be used as a proxy for the
simulation of STM experiments. The 3D LDOS file is then processed by a
specialized program **plstm** to produce a plot of the LDOS in
a 2D section at a given height in the unit cell (simulating the height
of a STM tip), or a simulated topography map by recording the z
coordinates with a given value of the LDOS.

The inputs to the STM workchain include a (possibly
already relaxed) structure and the protocol specification. The energy
window for the LDOS and the tip height or the LDOS iso-value can be in
principle specified by the user if full control is desired (probably
after evaluation of the results of the **SiestaBandsWorkchain**
workflow), but for the purposes of a turn-key solution, a range of
energies around the Fermi Level (or regions near to the HOMO and/or
LUMO), and a range of heights should automatically be selected by the
workflow and the results presented to the user for further
consideration. The workflow executes the **plstm** program via an
AiiDA plugin, which is also included in the **aiida-siesta**
distribution. Its parser stage returns an AiiDA 2D array whose
contents can be displayed by standard tools within AiiDA and the wider
Python ecosystem.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1 series, which
can be found in the development platform
(http://launchpad.net/siesta/).

Inputs
------

* **code**, a code associated to the Siesta plugin

* **stm_code**, a code associated to the STM (plstm)  plugin

* **structure**, class :py:class:`StructureData
  <aiida.orm.data.structure.StructureData>`

A structure. See the plugin documentation for more details.

* **protocol**, Str

Either "standard" or "fast" at this point.
Each has its own set of associated parameters.

- "standard"::
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.2,
                'dm_convergence_threshold': 1.0e-4,
                'forces_convergence_threshold': "0.02 eV/Ang",
                'min_meshcutoff': 100, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'md-type-of-run': "cg",
                'md-num-cg-steps': 10,
                'pseudo_familyname': 'lda-ag',
                # Future expansion. Add basis info, caveats, etc
                'atomic_heuristics': {
                    'H': { 'cutoff': 100 },
                    'Si': { 'cutoff': 100 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'DZP'
                }

- "fast"::
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.25,
                'dm_convergence_threshold': 1.0e-3,
                'forces_convergence_threshold': "0.2 eV/Ang",
                'min_meshcutoff': 80, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'md-type-of-run': "cg",
                'md-num-cg-steps': 8,
                'pseudo_familyname': 'lda-ag',
                # Future expansion. Add basis info, caveats, etc
                'atomic_heuristics': {
                    'H': { 'cutoff': 50 },
                    'Si': { 'cutoff': 50 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'SZP'
                }

The *atomic_heuristics* dictionary is intended to encode the
peculiarities of particular elements. It is work in progress.

The *basis* section applies globally for now.

Outputs
-------

* **output_structure** :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>` 

The final relaxed structure (if applicable)

* **stm_array** :py:class: `ArrayData <aiida.orm.data.array.ArrayData>` 

A 2D array holding the section or topography information.
  



