The protocols system
++++++++++++++++++++

Description
-----------

As explained in the corresponding section, the **SiestaBandsWorkchain** workflow 
can be used to perform simple DFT calculations, the relaxation of a structure and 
to visualize the electronic bands of the system under investigation. 
In order to use this tool, the user needs to manually select all the inputs of the workchain, 
being careful to pass the correct specifications to perform the task (see examples/workflos/example_first.py).
The package aiida_siesta provides also a set of pre-selected inputs to run a **SiestaBandsWorkchain**,
supporting the tasks of the relaxation of a structure and the calculations of bands.
In other words, the user can obtain a `builder` of the 
**SiestaBaseWorkChain** that is ready to be submitted. This `builder`, in fact, is pre-filled
with inputs selected according to very few options specified by the user.
The list of options is presented `here <how-to>`, however the main feature is the 
use of *protocols*. A *protocol* groups operational parameters for a Siesta calculation
and it is meant to offer a set of inputs with the desired balance of accuracy and efficiency.
At the moment only two protocols are shipped in the package, they are called 
*stringent_delta* and *standard_delta*. More on them is presented in the next to next subsection.
It is important to note that the implemented protocols are not, for the moment,
input parameters that are guaranteed to perform in any situation. They are only
based on reasonable assumptions and few tests. However, in the packag,e it is also implemented
a system that allows users to create their own protocols, as clarified in the next to next subsection.
Finally, it must be remembered that the `builder` produced according to a *protocol* and few other options is fully 
modifiable before submission, leaving full flexibility to the user.
We expect in the future to have more and more "know how" and improve the
reliability and richness of the available *protocols*.


Supported Siesta versions
-------------------------

The protocol system supports **only the MaX-1.0 version of Siesta**, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).


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

* *standard_delta*
  Pseudopotential ONCVPSPv0.4 (norm-conserving) of Pseudo Dojo in psml format, scalar relativistic,
  PBE and with *standard* accuracy. See Pseudo Dojo website for more info (link).
  Basis set apply globally, with size DZ and energy-shift of 100 meV. Meshcutoff is 100 Ry,
  electronic temp 25 meV, and a kpoint mesh with distance 0.2 are implemented.
  Concerning the trashold for convergence, we implement 1.e-3 tolerance for the density matrix,
  0.04 ev/ang for forces and 1 GPa for stress.
  This choice op parameters have been tested on crystal elements of the Delta test (link) and ... 

* *stringent_delta*
  Pseudopotential ONCVPSPv0.4 (norm-conserving) of Pseudo Dojo in psml format, scalar relativistic,
  PBE and with *stringent* accuracy. See Pseudo Dojo website for more info (link).
  Basis set apply globally, with size DZP and energy-shift of 100 meV. Meshcutoff is 500 Ry,
  electronic temp 25 meV, and a kpoint mesh with distance 0.062 are implemented.
  Concerning the trashold for convergence, we implement 1.e-4 tolerance for the density matrix,
  0.03 ev/ang for forces and 0.005 ev/ang^3 for stress.
  This choice op parameters have been tested on crystal elements of the Delta test (link) and ...

The management of the pseudos is, at the moment, very fragile. It imposes that the user
loads a pseudo_family with the correct name that is hard-coded for the each protocol.
This name are 'nc-sr-04_pbe_standard_psml' and 'nc-sr-04_pbe_stringent_psml' for *standard_delta* and
*stringent_delta* protocols respectively.
Therefore a user, before using protocol, needs to download the correct pseudos and
load them (verdi data psml uploadfamily) with the correct name.
---This last part will change soon, replaced with a proper setup-profile script ----

Few more variables are set for both plugins. They are related to mixing options: 
scf-mixer-history is set to 5, and scf-mixer-weight is 0.1. As only the Max-1.0 version 
of Siesta is supported, the default mixer is Pulay and the quantity mixed is the Hamiltonian.

The details explained above for the two implemented plugins refer to the global variables, however
each protocol implements also a dictionary called  *atomic_heuristics*. This dictionary is intended to encode the
peculiarities of particular elements. It is work in progress.

The full list of global parameters are collected in the `protocol_registry.yaml` file, located in 
aiida_siesta/workflows/utils. In there also the atom heuristics implemented can be explored.

How to create my protocols
--------------------------

The protocol system allows also to create customized protocol. In order to do that,
a file similar to `aiida_siesta/workflows/utils/protocol_registry.yaml`
must be created and here the custom protocols can be listed.
Then the path of this file must be added to the environment variable `AIIDA_SIESTA_PROTOCOLS`.
This will be sufficient to let aiida-siesta recognize the protocols.
Of course the file containing the new protocols must have the same structure of `protocol_registry.yaml`.
The protocol name should be the outer entry of the indentation.
Some keyword are mandatory, they are `description`, `parameters`, `basis` and `pseudo_family`. The `pseudo_family`
must contain the name of a family (Psml or Psf family) that has been already uploaded in the database.
The `parameters` and `basis` entries are transformed into dictionaries and directly passed
to AiiDA, but the syntax (lower case and '-' between words) must be respected as it is used
to introduce atom heuristic changes.
Optional entries are `kpoints` (distance and offset) and `relax` (containing the MD options,
it is very important to confine here the relaxation options and not mix them with other parameters,
this because there is a specific way to trigger relaxation in the protocol system, see next section).
If some atom heuristics are wanted, the user should know that the following are allowed:
for `parameters`, only the 'mesh-cutoff' is allowed (spin in the future) and the
parameter applies globally. This means that if we define a mesh-cutoff heuristic for element
"El", for a structure containing "El" the global 'mesh-cutoff' parameter will be substituted by
the heuristic one. This is meant to signal elements that requires a bigger 'mesh-cutoff' than normal.
If more elements in the structure have a custom 'mesh-cutoff', the biggest one will be selected.
For `basis`, we allow 'split-norm', 'polarization', 'energy-shift', 'size'. The 'split-norm' and 'energy-shift' apply
globally, exactly like 'mesh-cutoff'. The 'size' and' polarization' introduce instead a block
reporting the change of pao size and polarization schema only for the element under definition.

.. _how-to:

How to use protocols
--------------------

In this section we explain how to obtain a pre-filled builder according to a protocol,
that is ready to be submitted (or modified and then submitted).
Two are the ways:

* Make use of the class **BaseWorkChainInputsGenerator**.
  This class implements the method `get_builder` that return the pre filled `builder` for the 
  **SiestaBaseWorkChain**.
  Therefore something like::

        inp_gen = BaseWorkChainInputsGenerator()
        builder = inp_gen.get_builder(structure, calc_engines, protocol, path_generator, relaxation_type)
        #here user can modify builder befor submission.
        submit(builder)

  is sufficient to run the calculation. The parameters of `get_builder` are explained below.
  
* Call directly a filled builder from **SiestaBaseWorkChain**.
  Through the class `get_filled_builder`::

        builder = SiestaBaseWorkChain.get_filled_builder(structure, calc_engines, protocol)
        #here user can modify builder befor submission.
        submit(builder)

  The parameters of `get_filled_builder` of **SiestaBaseWorkChain** are the same of
  the `get_builder` of **BaseWorkChainInputsGenerator** and they are explained below.

In addition to protocol name, few more more parameters are necessary to `get_builder` of **BaseWorkChainInputsGenerator**
or `get_filled_builder` of **SiestaBaseWorkChain** in order to produce a correct builder. They are listed here:

* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`, *Mandatory*
  A structure. See the plugin documentation for more details.

* **calc_engine**, python `dict`, *Mandatory*
  A dictionary containing the specifications of the code to run. An example::

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

* **protocol**, python `str`, *Mandatory*
  The protocol name, selected among the available ones, as explained in the previous section.

* **path_generator**, python `str`. *Optional*
  The presence of this parameter triggers the calculation of bands.
  Two are the available value to pass as `path_generator`: "seekpath" or "legacy".
  They set the way the path in k-space is produced. This path is used to display the
  bands. While "seekpath" modify the structure running the calculation on an equivalent "conventional" 
  cell, "legacy" doesn't and preserves the input structure. However the "legacy" method is known to be wrong for 
  some particular structure.

* **relaxation_type**, python `str`. *Optional*
  The presence of this parameter triggers the possibility to relax the structure.
  The specifications of the relaxation_type are "atoms_only", "variable_cell" or "constant_volume",
  that should be self expalnatory.
  For the moment only the CG relaxation algorithm is implemented (in the future more will be added).

An example of the use is in aiida_siesta/examples/workflows/example_protocol.py.

We took the case of SiestaBaseWorkChain, but we implement a <WorkChain>InputsGenerator for each
<WorkChain> distributed in the package, this is highlighted in the corresponding subsection of the
<WorkChain> documentation.

Other methods of **BaseWorkChainInputsGenerator**
-------------------------------------------------

The only advantage of using the class **BaseWorkChainInputsGenerator** is that this class provides,
on top of `get_builder` explained above, a lot of methods that facilitate the task of exploring
the various options of the protocol system. For instance, there is a method listing all the available protocols,
the available relaxation types and so on.
