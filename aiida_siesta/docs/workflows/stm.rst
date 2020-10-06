STM workflow
++++++++++++

Description
-----------

The **SiestaSTMWorkchain** workflow consists in 3 steps:

* Performing of a siesta calculation on an input structure (including relaxation if needed) 
  through the **SiestaBaseWorkChain**.
* Performing of a further siesta calculation aimed to produce a .LDOS file.
* A call to the `plstm` code to post process the .LDOS file and
  create simulated STM images. The call is made via the
  **STMCalculation** plugin, which is also included in the ``aiida_siesta`` distribution.

The .LDOS file contains informations on the local density
of states (LDOS) in an energy window. The LDOS can be seen as a
"partial charge density" to which only those wavefunctions with
eigenvalues in a given energy interval contribute. In the
Tersoff-Hamann approximation, the LDOS can be used as a proxy for the
simulation of STM experiments. The 3D LDOS file is then processed by the
specialized program `plstm` to produce a 2D section in "constant-height" or 
"constant-current" mode, optionally projected on spin components
(see the header/manual for plstm, and note that non-collinear and spin-orbit 
modes are supported). 
The "constant-height" mode corresponds to the creation of 
a plot of the LDOS in a 2D section at a given height in the unit cell 
(simulating the height of a STM tip). The "constant-current" mode
simulates the topography map by recording the z
coordinates with a given value of the LDOS.

The inputs to the STM workchain include all the inputs of the **SiestaBaseWorkChain**
to give full flexibility on the choice of the siesta calculation
parameters. The energy window for the LDOS is specified respect to the Fermi energy.
In fact, a range of
energies around the Fermi Level (or regions near to the HOMO and/or
LUMO) are the meaninful energies for the STM images production. 
The tip height ("constant-height" mode) or the LDOS iso-value ("constant-current" mode)
must be specified by the user in input.
The workchain returns an AiiDA ArrayData object whose
contents can be displayed by standard tools within AiiDA and the wider
Python ecosystem.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the      
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.



Inputs
------

* All the inputs of the **SiestaBaseWorkChain**, as explained
  :ref:`here <siesta-base-wc-inputs>`.

.. |br| raw:: html

    <br />

* **stm_code**, class :py:class:`Code  <aiida.orm.Code>`, *Mandatory*

  A code associated to the STM (plstm) plugin (siesta.stm). See plugin documantation for more details.

.. |br| raw:: html

    <br />

* **stm_mode**, class :py:class:`Str <aiida.orm.Str>`, *Mandatory*

  Allowed values are ``constant-height`` or ``constant-current``, corresponding to the two
  operation modes of the STM that are supported by the plstm code.

.. |br| raw:: html

    <br />


* **stm_value**, class :py:class:`Float <aiida.orm.Float>`, *Mandatory*

  The value of height or current at which the user wants to simulate the
  STM. This value represents the tip height in "constant-height" mode
  or the LDOS iso-value in "constant-current" mode.
  The height must be expressed in `Angstrom`, the current in `e/bohr**3`.

.. |br| raw:: html

    <br />


* **emin**, class :py:class:`Float  <aiida.orm.Float>`, *Mandatory*

  The lower limit of the energy window for which the LDOS is to be
  computed (in eV and respect to the Fermi level).

.. |br| raw:: html

    <br />

* **emax**, class :py:class:`Float <aiida.orm.Float>`, *Mandatory*

  The upper limit of the energy window for which the LDOS is to be
  computed (in `eV` and respect to the Fermi level).

.. |br| raw:: html

    <br />

* **stm_spin**, class :py:class:`Str <aiida.orm.Str>`, *Mandatory*

  Allowed values are ``none``, ``collinear`` or ``non-collinear``.
  Please note that this keyword only influences the STM post process!
  It does not change the parameters of the siesta calculation, that must
  be specified in the **parameters** input port.
  In fact, this keyword will be automatically reset if a `stm_spin`
  option incompatible with the parent siesta spin option is chosen.
  A warning will be issued in case this happens.
  This keyword also influences the structure of the output port
  **stm_array**. If fact, if the ``non-collinear`` value is chosen, the
  workflow automatically performs the STM analysis in the three
  spin components and for the total charge option, resulting in a
  richer **stm_array** (see description in the Outputs section).

.. |br| raw:: html

    <br />

* **stm_options**, class :py:class:`Dict <aiida.orm.Dict>`, *Optional*
  
  This dictionary can be used to specify the computational resources to
  be used for the STM calculation (the `plstm` code). It is optional
  because, if not specified, the same resources of the siesta calculations
  are used, except that the parallel options are stripped off.
  In other words, by default, the `plstm` code runs on a single processor. 


Outputs
-------

* **stm_array** :py:class:`ArrayData <aiida.orm.ArrayData>` 

  In case the **stm_spin** is ``none`` or ``collinear`` this output port
  is a collection of three 2D arrays (`grid_X`, `grid_Y`, `STM`) holding the section or
  topography information. Exactly like the output of the STM plugin.
  In case the **stm_spin** is ``non-collinear``, this output port
  is a collection of six 2D arrays (`grid_X`, `grid_Y`, `STM_q`, `STM_sx`, `STM_sy`, `STM_sz`)
  holding the section or topography information for the total charge STM analysis and 
  the three spin components.
  Both cases follow the `meshgrid` convention in
  Numpy. A contour plot can be generated with the `get_stm_image.py`
  script in the repository of examples. The `get_stm_image.py` script
  automatically detects how many arrays are in **stm_array**, therefore it is 
  completely general.

.. |br| raw:: html

    <br />

* **output_structure** :py:class:`StructureData <aiida.orm.StructureData>`

  Present only if the siesta calculation is moving the ions.  Cell and ionic
  positions refer to the last configuration, on which the STM analysis is performed.

  

Protocol system
---------------

The protocol system is available for this WorkChain. The ``SiestaSTMWorkchain.inputs_generator()``
makes available all the methods explained in the :ref:`protocols documentation <how-to>`, but
``get_filled_builder`` now requires in inputs also the ``stm_mode`` (a python `str <str>`, accepted values 
are "constant-height" and "constant-current") and ``stm_value`` (a python `float <float>` indicating
the value of height in Ang or current in e/bohr**3).
Moreover in the ``calc_engines`` dictionary, also indications on the resources for the stm calculation must
specified, following the syntax of this example::

   calc_engines = {
     'siesta': {
         'code': codename,
         'options': {'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1}, "max_wallclock_seconds": 3600 }
         },
     'stm': {
         'code': stmcodename,
         'options': {'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1}, "max_wallclock_seconds": 1360 }
         }
     }

The STM spin mode is chosen accordingly to the ``spin`` input passed to ``get_filled_builder``,
setting "collinear" stm_spin in case of polarized calculation, "non-collinear" in case of 
"spin-orbit" or "non-collinear" calculations and no spin in case of an unpolarized calculation.
Therefore, if, for instance, the user wants to post-process a spin calculation with "no-spin"
STM mode, he/she needs to manually modify the builder before submission.
Also the **emin** and **emax** inputs of **SiestaSTMWorkchain** are internally chosen
by the inputs generator: they select an energy window of `6 eV` below the Fermi energy.
If the choice doesn't suit the purpose, the user can manually modify the builder before
submission.
