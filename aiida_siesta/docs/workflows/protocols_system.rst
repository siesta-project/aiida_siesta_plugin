The protocols system
++++++++++++++++++++

Description
-----------

As explained in the corresponding sections, the **SiestaBandsWorkchain** workflow
is the suggested tool to run siesta calculations in aiida.
In order to use this tool, the user needs to manually select all the inputs of the workchain, 
being careful to pass the correct specifications to perform the calculation
(see examples/workflos/example_first.py).
The package aiida_siesta provides also a set of pre-selected inputs to run a **SiestaBandsWorkchain**,
and the other WorkChains,
supporting the tasks of the relaxation of a structure and the calculations of bands.
In other words, the user can obtain a `builder` of the 
**SiestaBaseWorkChain** that is ready to be submitted. This `builder`, in fact, is pre-filled
with inputs selected according to the structure under investigation and very few options specified by the user.
The list of options is presented :ref:`here <how-to>`, however the main feature is the 
use of *protocols*. A *protocol* groups operational parameters for a Siesta calculation
and it is meant to offer a set of inputs with the desired balance of accuracy and efficiency.
At the moment only two protocols are shipped in the package, they are called 
*stringent* and *standard*. More on them is presented in the next to next subsection.
It is important to note that the implemented protocols are not, for the moment,
input parameters that are guaranteed to perform in any situation. They are only
based on reasonable assumptions and few tests. However, in the package it is also implemented
a system that allows users to create their own protocols, as clarified :ref:`here <custom-prot>`.
Finally, it must be remembered that the `builder` produced according to a *protocol* and few other options is fully 
modifiable before submission, leaving full flexibility to the user.
We expect in the future to have more and more "know how" and improve the
reliability and richness of the available *protocols*.
We focus here on the description of the use of protocols for the **SiestaBandsWorkchain**,
but the same system is available for all the workchains distributed in this package.
A small paragraph in the documentation of each WorkChain will explain the details of
the usage of protocols for that particular workchain.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).

.. The protocol system supports **only the MaX-1.0 version of Siesta**, which
.. can be found in the development platform
.. (https://gitlab.com/siesta-project/siesta).


Available protocols
-------------------

With the word *protocol* we mean a series of suggested inputs for AiiDA
WorkChains that allow users to more easily automatize their workflows.
These inputs reflects a certain set of operational parameters for a Siesta
calculation. The choice of the inputs of a DFT simulation should be carefully tested
for any new system. Therefore the use of protocols, in place of a careful and tested
choice of inputs, it is always somehow limiting. It can be, however, 
considered a good starting point.
This is the very beginning of the development and, for the moment, only
two very basic protocols are implemented.
A description of the global variables of each protocol are now reported.

* *standard*

  Pseudopotential ONCVPSPv0.4 (norm-conserving) of Pseudo Dojo in psf format, scalar relativistic,
  PBE and with *standard* accuracy. Download at https://icmab.es/leem/SIESTA_MATERIAL/tmp_PseudoDojo/nc-sr-04_pbe_standard-psf.tgz.
  Basis set apply globally, with size DZ and energy-shift of 100 meV. Meshcutoff is 100 Ry,
  electronic temp 25 meV, and a kpoint mesh with distance 0.2 are implemented.
  Concerning the trashold for convergence, we implement 1.e-3 tolerance for the density matrix,
  0.04 ev/ang for forces and 1 GPa for stress.
  This choice of inputs (plus some atom heuristics - see below) have been run for a all
  the crystal elements up to the element Po (excluding lanthanides) but performances have not been tested.

.. |br| raw:: html

    <br />


* *stringent*

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
This name are 'nc-sr-04_pbe_standard-psf' and 'nc-sr-04_pbe_stringent-psf' for *standard* and
*stringent* protocols respectively.
Therefore a user, before using protocol, needs to download the correct pseudos and
load them (verdi data psf uploadfamily) with the correct name.
---This last part will change soon, replaced with a proper setup-profile script ----

Few more variables are set for both protocols. They are related to mixing options: 
scf-mixer-history is set to 5, and scf-mixer-weight is 0.1. As only the Max-1.0 version 
of Siesta is supported, the default mixer is Pulay and the quantity mixed is the Hamiltonian.

The details explained above for the two implemented plugins refer to the global variables, however
each protocol implements also a dictionary called *atomic_heuristics*. This dictionary is intended to encode the
peculiarities of particular elements. It is work in progress.

The full list of global parameters are collected in the `protocol_registry.yaml` file, located in 
aiida_siesta/workflows/utils. In there also the atom heuristics implemented can be explored.


.. _how-to:

How to use protocols
--------------------

In this section we explain how to obtain a pre-filled builder according to a protocol
and an input structure, that is ready to be submitted (or modified and then submitted).

First of all, the 'nc-sr-04_pbe_standard-psf' and 'nc-sr-04_pbe_stringent-psf' set of
pseudopotentials must be downloaded and stored in the database in a family
with the same name::
        
        wget https://icmab.es/leem/SIESTA_MATERIAL/tmp_PseudoDojo/nc-sr-04_pbe_standard-psf.tgz
        wget https://icmab.es/leem/SIESTA_MATERIAL/tmp_PseudoDojo/nc-sr-04_pbe_stingent-psf.tgz
        tar -xf nc-fr-04_pbe_standard-psf.tgz
        tar -xf nc-sr-04_pbe_stringent-psf.tgz
        verdi data psf uploadfamily nc-fr-04_pbe_standard-psf nc-fr-04_pbe_standard-psf "Scalar-relativistic psf standard"
        verdi data psf uploadfamily nc-sr-04_pbe_standard-psf nc-sr-04_pbe_stringent-psf "Scalar-relativistic psf stringent"


Once this first step is done,
the pre-filled builder can be
accessed through the method `inputs_generator` of the **SiestaBaseWorkChain**, like 
in this example::

        inp_gen = SiestaBaseWorkChain.inputs_generator()
        builder = inp_gen.get_filled_builder(structure, calc_engines, protocol)
        #here user can modify builder befor submission.
        submit(builder)

The parameters of `get_filled_builder` of **SiestaBaseWorkChain** are explained here:

* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`, *Mandatory*

  A structure. See the plugin documentation for more details.

.. |br| raw:: html

    <br />

* **calc_engine**, python `dict`, *Mandatory*

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

  The dictionary must present `siesta` as upper level key of the dictionary. This might seem unnecessary, but
  will become fundamental for the use of protocols in more complicated workchain, involving not only
  the siesta plugin, but also, for instance, the stm plugin.

.. |br| raw:: html

    <br />

* **protocol**, python `str`, *Mandatory*

  The protocol name, selected among the available ones, as explained in the previous section.

.. |br| raw:: html

    <br />

* **bands_path_generator**, python `str`, *Optional*

  The presence of this parameter triggers the calculation of bands.
  Two are the available value to pass as `bands_path_generator`: "seekpath" or "legacy".
  They set the way the path in k-space is produced. This path is used to display the
  bands. While "seekpath" modify the structure running the calculation on an equivalent "conventional" 
  cell, "legacy" doesn't and preserves the input structure. However the "legacy" method is known to 
  have bugs for certain structure cells.

.. |br| raw:: html

    <br />

* **relaxation_type**, python `str`, *Optional*

  The presence of this parameter triggers the possibility to relax the structure.
  The specifications of the relaxation_type are "atoms_only", "variable_cell" or "constant_volume",
  that should be self expalnatory.
  For the moment only the CG relaxation algorithm is implemented (in the future more will be added).

.. |br| raw:: html

    <br />

* **spin**, python `str`, *Optional*

  The presence of this parameter triggers the spin options.
  The specifications of the spin are the one of modern version of Siesta, they are
  "polarized", "non-collinear" and "spin-orbit".
  If no spin option is defined, the calculation will not be spin polarized.

An example of the use is in aiida_siesta/examples/workflows/example_protocol.py.

The method `get_filled_builder` is definitely the most important tool offered by the `inputs_generator`,
however through this property of **SiestaBaseWorkChain** other methods that facilitate the task of exploring
the various options of the protocol system are available. For instance, there is a method listing all the available protocols,
the available relaxation types and so on.

.. _custom-prot:

How to create my protocols
--------------------------

The protocol system allows also to create customized protocol. To this end, a
file similar to `aiida_siesta/workflows/utils/protocol_registry.yaml`
must be created, listing the custom protocols.
Then the path of this file must be added to the environment variable `AIIDA_SIESTA_PROTOCOLS`.
This will be sufficient to let aiida-siesta recognize the protocols.

The file containing the customized protocols must have the same structure of `protocol_registry.yaml`.
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
texactly for this purpose.
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
For `basis`, we allow 'split-tail-norm', 'polarization' and 'size'. The 'size' and' polarization' introduce a block
reporting the change of pao size and polarization schema only for the element under definition.
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
