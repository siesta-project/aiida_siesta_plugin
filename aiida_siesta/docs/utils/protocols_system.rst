The protocols system
++++++++++++++++++++

Description
-----------

In order to submit **SiestaCalculations**, the user needs to manually select all the inputs, 
being careful to pass the correct specifications to perform the calculation
(as explained in the :ref:`corresponding section <siesta-plugin-inputs>`).
The package ``aiida_siesta`` provides also a set of pre-selected inputs to run a **SiestaCalculations**,
and the WorkChains distributed in the package,
supporting the tasks of the relaxation of a structure and the calculations of bands.
In other words, the user can obtain a ``builder`` of the 
**SiestaCalculation** that is ready to be submitted. This ``builder``, in fact, is pre-filled
with inputs selected according to the structure under investigation and very few options specified by the user.
The lengthy inputs selection is substitute by::

        inp_gen = SiestaCalculation.inputs_generator()
        builder = inp_gen.get_filled_builder(structure, calc_engines, protocol)


The list of options to obtain the builder is presented :ref:`here <how-to>`, however the main feature is the 
use of *protocols*. A *protocol* groups operational parameters for a Siesta calculation
and it is meant to offer a set of inputs with the desired balance of accuracy and efficiency.
At the moment only one protocol is shipped in the package, it is called 
*standard_psml*. More on it is presented in the next to next subsection.
It is important to note that the implemented protocols are not, for the moment,
input parameters that are guaranteed to perform in any situation. They are only
based on reasonable assumptions and few tests. However, in the package it is also implemented
a system that allows users to create their own protocols, as clarified :ref:`here <custom-prot>`.
Finally, it must be remembered that the ``builder`` produced according to a *protocol* and few other options is fully 
modifiable before submission, leaving full flexibility to the user.
We expect in the future to have more and more "know how" and improve the
reliability and richness of the available *protocols*.

We focus here on the description of the use of protocols for the **SiestaCalculation**,
but the same system is available for all the WorkChains distributed in this package.
A small paragraph in the documentation of each WorkChain will explain the details of
the usage of protocols for that particular WorkChain.


Supported Siesta versions
-------------------------

.. At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
.. can be found in the development platform
.. (https://gitlab.com/siesta-project/siesta).

The protocol system, at the moment, requires a version of siesta 
with support for psml pseudopotential. At least **the MaX-1.0 release of Siesta**, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta), meets this requirement.
For more up to date info on compatibility, please check the
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.



Available protocols
-------------------

With the word *protocol* we mean a series of suggested inputs for AiiDA
CalJobs/WorkChains that allow users to more easily automatize their workflows.
These inputs reflects a certain set of operational parameters for a Siesta
calculation. The choice of the inputs of a DFT simulation should be carefully tested
for any new system. Therefore the use of protocols, in place of a careful and tested
choice of inputs, it is always somehow limiting. It can be, however, 
considered a good starting point.
This is the very beginning of the development and, for the moment, only
one very basic protocol is implemented.
A description of its variables is now reported. Each protocol contain a section
with global variables and an *atomic_heuristics* dictionary, a dictionary intended to encode the
peculiarities of particular elements.


* *standard_psml*

  The full list of variables for this protocol are collected in the `protocol_registry.yaml` file, located in
  ``aiida_siesta/utils``.

  * *global variables*

    Pseudopotential ONCVPSPv0.4 (norm-conserving) of Pseudo Dojo in psml format, scalar relativistic,
    PBE and with *standard* accuracy (download available from the `PseudoDojo`_ web site).
    Basis set apply globally, with size DZP and ``energy-shift`` of 50 meV. The ``mesh-cutoff`` is 200 Ry,
    ``electronic-temp`` 25 meV, and a kpoint mesh with distance 0.1 are implemented.
    Concerning the trashold for convergence, we implement 1.e-4 tolerance for the density matrix,
    0.04 ev/ang for forces and 0.1 GPa for stress.
    Few more global variables are related to mixing options:
    ``scf-mixer-history`` is set to 5, and ``scf-mixer-weight`` is 0.1. As only the Max-1.0 version
    of Siesta is supported, the default mixer is Pulay and the quantity mixed is the Hamiltonian.

  .. |br| raw:: html

    <br />

 
  * *atomic_heuristics*

    The element "Ag" requires a bigger ``mesh-cutoff`` because ``mesh-cutoff = 200 Ry`` was leading to a
    "Failure to converge standard eigenproblem" error for the "Ag" elemental crystal.
    Custom basis for "Ca","Sr","Ba" are necessary because the automatic generation results
    in a too-large radius for the "s" orbitals. The "Hg" custom basis introduces an increment of
    all radii of 5% compared to the automatic generated orbitals and adds a Z orbital for the "p"
    channel, while removing polarization.
    The elements "Li", "Be", "Mg", Na", "Fe", "Mn", "Sb" require a bigger 
    ``mesh-cutoff`` because ``mesh-cutoff = 200 Ry`` resulted in
    a discontinuous equation of state.

  This choice of parameters have been tested on crystal elements up to the 
  element "Rn" and compared with the reference equation of state of the
  `DeltaTest`_ project, resulting on an average delta value of 7.1 meV.
  The parameters of this protocol for noble gasses do not result in an a minimum of the equation of state.
  Because Van der Waals forces are not included in the calculation, the result is not surprising.
  We warn users to use with care this protocol for noble gasses.
  It is important to stress that the present protocol has not been conceived to produce
  good results for the Delta test; the basis sets are mostly automatic and the choice of
  mesh-cutoff / kpoints-mesh is farely loose. The average value for the delta (7.1 meV)
  is just an indication that the parameters' choice gives reasonable results for elemental crystals.
  We are working on a more accurate (and expensive) protocol that will provide much better
  values of delta.
  New tests and checks on the *standard_psml* protocol will be added in the aiida-siesta 
  `wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Protocols-validations>`_.


  
.. Maximum delta is 28 meV for "Ne" and "Ar".
  
..  Download at https://icmab.es/leem/SIESTA_MATERIAL/tmp_PseudoDojo/nc-sr-04_pbe_standard-psf.tgz.
  Basis set apply globally, with size DZ and energy-shift of 100 meV. Meshcutoff is 100 Ry,
  electronic temp 25 meV, and a kpoint mesh with distance 0.2 are implemented.
  Concerning the trashold for convergence, we implement 1.e-3 tolerance for the density matrix,
  0.04 ev/ang for forces and 1 GPa for stress.
  This choice of inputs (plus some atom heuristics - see below) have been run for a all
  the crystal elements up to the element Po (excluding lanthanides) but performances have not been tested.
  Pseudopotential ONCVPSPv0.4 (norm-conserving) of Pseudo Dojo in psf format, scalar relativistic,
  PBE and with *stringent* accuracy. Download at https://icmab.es/leem/SIESTA_MATERIAL/tmp_PseudoDojo/nc-sr-04_pbe_standard-psf.tgz.
  Basis set apply globally, with size DZP and energy-shift of 50 meV. Meshcutoff is 500 Ry,
  electronic temp 25 meV, and a kpoint mesh with distance 0.062 are implemented.
  Concerning the trashold for convergence, we implement 1.e-4 tolerance for the density matrix,
  0.01 ev/ang for forces and 0.05 GPa for stress.
  This choice of parameters (plus some atom heuristics - see below)
  have been tested on crystal elements up to the element Au (excluding
  lanthanides and noble gasses) and compared with the reference equation of state of the
  `DeltaTest`_ project, resulting in values of delta below 10 meV for all elements except
  "N", "Ca", "Ga", "Ge", "As", "Sr", "In", "Sb", "Ba".
  Investigations are on-going in order to improve the performance of the available set, however
  it must be remembered that the test on crystal elements has very limited meaning
  when the atoms are in more complex chemical environments.

The management of the pseudos is, at the moment, very fragile. It imposes that the user
loads a pseudo_family with the correct name that is hard-coded for the each protocol.
This name is 'nc-sr-04_pbe_standard_psml' for the *standard_psml* protocol.
Therefore a user, before using protocol, needs to download the correct pseudos and
load them (see next section) with the correct name.
---This last part will change soon, replaced with a proper setup-profile script ----

.. _how-to:

How to use protocols
--------------------

In this section we explain how to obtain a pre-filled builder according to a protocol
and an input structure, that is ready to be submitted (or modified and then submitted).

First of all, the 'nc-sr-04_pbe_standard_psml' set of
pseudopotentials must be downloaded from `PseudoDojo`_ and stored in the database in a family
with the same name. From command line::
     
      wget http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard_psml.tgz
      mkdir nc-sr-04_pbe_standard_psml
      tar -xf nc-sr-04_pbe_standard_psml.tgz -C nc-sr-04_pbe_standard_psml
      verdi data psml uploadfamily nc-sr-04_pbe_standard_psml nc-sr-04_pbe_standard_psml "Scalar-relativistic psf standard"
        
..      wget https://icmab.es/leem/SIESTA_MATERIAL/tmp_PseudoDojo/nc-sr-04_pbe_standard-psf.tgz
        tar -xf nc-sr-04_pbe_standard-psf.tgz
        verdi data psf uploadfamily nc-sr-04_pbe_standard-psf nc-sr-04_pbe_stringent-psf "Scalar-relativistic psf stringent"


Once this first step is done, the pre-filled builder can be
accessed through the method ``inputs_generator`` of the **SiestaCalculation**
(and of any other workchain). 
For example::

        from aiida_siesta.calculations.siesta import SiestaCalculation
        inp_gen = SiestaCalculation.inputs_generator()
        builder = inp_gen.get_filled_builder(structure, calc_engines, protocol)
        #here user can modify builder befor submission.
        submit(builder)

The arguments of ``get_filled_builder`` of the input generator are explained here:

* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`, *Mandatory*

  A structure. See the :ref:`plugin documentation <siesta-plugin-inputs>` for more details.

.. |br| raw:: html

    <br />

* **calc_engine**, python :py:class:`dict`, *Mandatory*

  A dictionary containing the specifications of the code to run and the computational
  resources. An example::

        calc_engines = {
            'siesta': {
                'code': codename,
                'options': {
                        'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1},
                        'max_wallclock_seconds': 360, 
                        'queue_name': 'DevQ', 
                        'withmpi': True, 
                        'account': "tcphy113c"
                 }
            }
        }

  The dictionary must present ``siesta`` as upper level key of the dictionary. This might seem unnecessary, but
  will become fundamental for the use of protocols in more complicated WorkChain, involving not only
  the siesta plugin, but also, for instance, the stm plugin. The structure of ``calc_engines`` for each
  WorkChain of the package will be specified in the WorkChain documentation.

.. |br| raw:: html

    <br />

* **protocol**, python :py:class:`str`, *Mandatory*

  The protocol name, selected among the available ones, as explained in the previous section.

.. |br| raw:: html

    <br />

* **bands_path_generator**, python :py:class:`str`, *Optional*

  The presence of this parameter triggers the calculation of bands.
  Two are the available value to pass as `bands_path_generator`: "seekpath" or "legacy".
  They set the way the path in k-space is produced. This path is used to display the
  bands. While "seekpath" modify the structure running the calculation on an equivalent "conventional" 
  cell, "legacy" doesn't and preserves the input structure. However the "legacy" method is known to 
  have bugs for certain structure cells.

.. |br| raw:: html

    <br />

* **relaxation_type**, python :py:class:`str`, *Optional*

  The presence of this parameter triggers the possibility to relax the structure.
  The specifications of the relaxation_type are "atoms_only", "variable_cell" or "constant_volume",
  that should be self expalnatory.
  For the moment only the CG relaxation algorithm is implemented (in the future more will be added).

.. |br| raw:: html

    <br />

* **spin**, python :py:class:`str`, *Optional*

  The presence of this parameter triggers the spin options.
  The specifications of the spin are the one of modern version of Siesta, they are
  "polarized", "non-collinear" and "spin-orbit".
  If no spin option is defined, the calculation will not be spin polarized.

An example of the use is in `aiida_siesta/examples/plugins/siesta/example_protocol.py`

The method ``get_filled_builder`` is definitely the most important tool offered by the ``inputs_generator``,
however through the ``inputs_generator`` other methods can be accessed to explore
the various options of the protocol system. For instance, there is a method listing all the available protocols,
the available relaxation types and so on.

.. _custom-prot:

How to create my protocols
--------------------------

The protocol system allows also to create customized protocol. To this end, a
file similar to `aiida_siesta/utils/protocol_registry.yaml`
must be created, listing the custom protocols.
Then the path of this file must be added to the environment variable `AIIDA_SIESTA_PROTOCOLS`.
This will be sufficient to let aiida-siesta recognize the protocols.
The file containing the customized protocols must have the same structure of `protocol_registry.yaml`.
An example::

        my_protocol:
          description: 'My description'
          parameters:
            xc-functional: "GGA"
            xc-authors: "PBE"
            mesh-cutoff: '200 Ry'
            ...
          spin_additions:
            write-mulliken-pop: 1
          relax_additions:
            scf-dm-tolerance: 1.e-4
            md-max-force-tol: '0.04 eV/ang'
            md-max-stress-tol: '0.1 GPa'
          basis:
            pao-energy-shift: '50 meV'
            pao-basis-size: 'DZP'
          pseudo_family: 'nc-sr-04_pbe_standard_psml'
          kpoints:
            distance: 0.1
            offset: [0., 0., 0.]
          atomic_heuristics:
            Li:
              parameters:
                mesh-cutoff: '250 Ry'
              basis:
                polarization: 'non-perturbative'
                pao-block: "Li 3 \n  ... "
                split-tail-norm: True

The protocol name should be the outer entry of the indentation.
For each protocol, some keyword are mandatory. They are `description`, `parameters`, `basis` and `pseudo_family`. 
The `pseudo_family`
must contain the name of a family (Psml or Psf family) that has been already uploaded in the database.
The number of elements covered by your pseudo family will limit the materials you
can simulate with your protocol.
The `parameters` and `basis` entries are transformed into dictionaries and passed
to AiiDA after possible modifications due to atom heuristics or spin/relax additions.
For this reason, the syntax (lower case and '-' between words) must be respected in full.

Two optional keywords are `relax_additions` and `spin_additions`.
This two entries are not meant to host the siesta keywords that activate the relaxation or spin options,
but possible additions/modifications to the `parameters` entry, to apply in case of relaxation
or spin. When the use of protocols is called and the relax/spin options are requested (see `here <how-to>`_),
the system will automatically take care of introducing the correct siesta keyword (`MD.TypeOfRun`, 
`MD.VariableCell`, `spin` etc.) that are indispensable to run the task. However, it might happen that
a user desires a more loose `scf-dm-tolerance` for the task of the relaxation or a different `scf-mixer-weight`
when the spin is active. The `relax_additions` and `spin_additions` keywords have been created
exactly for this purpose.
Please be carefull that (except for the `mesh-cutoff`) if a keyword in `spin_additions` or 
`relax_additions` is already present in `parameters`, its value in `parameters` will overriden.
In other words, values in `spin_additions` or `relax_additions` have priority compared to the one
in `parameters`. Moreover `relax_additions` has priority respect to `spin_additions`.
For the `mesh-cutoff` the situation is different, because the biggest value will always be
considered, no metter where it is specified.
Another optional entry is `kpoints`, where a `distance` and an `offset` only can be specified.
The system will take care to create a uniform mesh for the structure under investigation with
a density that correspond to a distance (in 1/Angstrom) between adjacent kpoints equal to `dinstance`.

The final allowed (optional) keyword is `atomic_heuristics`. 
In it, two only sub-keys are allowed: `parameters` and `basis`.
In `parameters`,  only a 'mesh-cutoff' can be specified. This `mesh-cutoff` applies globally
and only if it is the biggest one among the all `mesh-cutoff` that apply.
This system is meant to signal elements that requires a bigger 'mesh-cutoff' than normal.
For `basis`, we allow 'split-tail-norm', 'polarization', 'size' and 'pao-block'. The 'size' and' polarization' introduce a block
reporting the change of pao size and polarization schema only for the element under definition.
The 'pao-block' allows to specify an explicit "block Pao-basis" for the element.
The 'split-tail-norm' instead activate in siesta the key 'pao-split-tail-norm', that applies globally.

We conclude this subsection with few more notes to keep in mind. First, the units mut be specified for each siesta keyword
that require units and they must be consisten throughout the protocol. This means that it is not possible
to define 'mesh-cutoff' in Ry in `parameters`, but in eV in the `atomic_heuristics`.
Second, it is up to the creator to remember to introcude the correct 'xc-functional' and 'xc-authors'
keywords in the protocol that matches the same exchange-correlation functional of the pseudos in the
pseudo family. This also means that we do not support pseudos presenting
different exchange-correlation functionals in the same family. Third, we impose a description for
each protocol because in the description the creator must underline the limitations of the protocol.
For instance, the case when a certain protocol do not support spin-orbit as the pseudos are not relativistics.
The schema we presented here is certanly not perfect and it is far to cover all the possible situations,
however it must be remembered that any user has always the chance to modify the inputs (builder) before submission.

.. _DeltaTest: https://molmod.ugent.be/deltacodesdft
.. _PseudoDojo: http://www.pseudo-dojo.org/
