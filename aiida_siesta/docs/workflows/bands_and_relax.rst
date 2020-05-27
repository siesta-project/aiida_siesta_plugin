Bands and relaxation with protocols
+++++++++++++++++++++++++++++++++++

Description
-----------

As explained in the corresponding section, the **SiestaBandsWorkchain** workflow 
can be used to perform the relaxation of a structure and to visualize the
electronic band structure of the system under investigation. In that case, the user
needs to manually select all the inputs of the workchain, being careful
to pass the correct specifications to perform the task.
The package aiida_siesta provides also a set of preselected inputs
for the tasks of the relaxation of a structure and the calculations of bands.
This is implemented through two classes, the class **SiestaRelaxInputsGenerator**
and **SiestaBandsInputsGenerator**. Both this classes implement a method `get_builder`
that return a `builder` for the **SiestaBaseWorkChain** with pre-compiled inputs,
selected according to some options passed to `get_builder`. The available options and
their specifications are explained below, separately for the task of relaxation and of 
the bands. However one common feature is implemented for both cases: the 
use of *protocols*. A *protocol* groups operational parameters for a Siesta calculation
and it is meant to offer a set of inputs with the desired balance of accuracy and efficiency.
At the moment only two protocols are implemented, they are called 
*stringent_delta* and *standard_delta*. More on them is presented in the next subsection.
It is important to note that the implemented protocols are not, for the moment,
tested input parameters that are guaranteed to perform in any situation. They are only
based on reasonable assumptions. However it must be remembered that the 
`builder` produced according to a *protocol* and few other options is fully 
modifiable before submission, leaving full flexibility to the user.
We expect in the future to have more and more "know how" and improve the
reliability and richness of the available *protocols*.


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
A description of the global variables of each protocol are reported below, however 
each protocol implements also a dictionary called  *atomic_heuristics*. This dictionary is intended to encode the
peculiarities of particular elements. It is work in progress.

* *standard_delta*
  Pseudopotential ONCVPSPv0.4 (norm-conserving) of Pseudo Dojo in psml format, scalar relativistic,
  PBE and with *standard* accuracy. See Pseudo Dojo website for more info (link).
  Basis set apply globally, with size DZ and energy-shift of 100 meV. Meshcutoff is 100 Ry,
  electronic temp 25 meV, and a kpoint mesh with dinstance 0.2 are implemented.
  Concerning the trashold for convergence, we implement 1.e-3 tolerance for the density matrix,
  0.04 ev/ang for forces and 1 GPa for stress.
  This choice op parameters have been tested on crystal elements of the Delta test (link) and ... 

* *stringent_delta*
  Pseudopotential ONCVPSPv0.4 (norm-conserving) of Pseudo Dojo in psml format, scalar relativistic,
  PBE and with *stringent* accuracy. See Pseudo Dojo website for more info (link).
  Basis set apply globally, with size DZP and energy-shift of 100 meV. Meshcutoff is 500 Ry,
  electronic temp 25 meV, and a kpoint mesh with dinstance 0.062 are implemented.
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

Few more variable set for both plugins, pulay options: num pulay 4, mixing 0.3. The rest is siesta default

Relaxation
----------

In addition to protocols, few more inputs are necessary to `get_builder` of **SiestaRelaxInputsGenerator**
in order to produce a correct builder:

* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`
  A structure. See the plugin documentation for more details.

* **calc_engine**, python dict
  A dictionary containing the specifications of the code to run. An example::

        calc_engines = {
            'relaxation': {
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

  The dictionary must present the key `relaxation` as upper level dictionary. This is not strictly necessary, but 
  it is used in this way in order to show that such schema could be extended to task where more than one
  step is necessary and a different code is necessary for each step.

* **relaxation_type**, python str.
  The specification of the relaxation_typ, meaning if "atom_only", "variable_cell" or "constant_volume"
  For the moment only CG algorith is implemented (in the future more will be added)

* **protocol**, python str.
  Explained above

* **threshold_forces**, python float. *Optional*
  If set, in ev/ang, overried the default of the protocol.

* **threshold_stress**, python float. *Optional*
  If set, in ev/ang^3, overried the default of the protocol


Therefore it is sufficient to have something like::

        rel_inp_gen = SiestaRelaxInputsGenerator()
        builder = rel_inp_gen.get_builder(
                structure=structure, calc_engines=calc_engines, protocol=protocol, 
                relaxation_type=relaxation_type)

To obtain a builder that can be directly submitted (`submit(builder)`) and it performs the task of relaxing a structure.
For the outputs, see specifications of the **SiestaBandsWorkchain**. Any input can be modified manually before
submission as the builder is not a stored object.
An example of the use is in aiida_siesta/examples/workflows/relax_with_protocol.py.

Bands
-----

In addition to protocols, few more inputs needs to be passed to `get_builder`:

* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`
  A structure. See the plugin documentation for more details.

* **calc_engine**, python dict
  A dictionary containing the specifications of the code to run. An example::

        calc_engines = {
            'bands': {
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

  The dictionary must present the key `bands` as upper level dictionary.

* **path_generator**, python str.
  The specification of the way the path in k-space is produced. This path is used to display the
  bands. Two options are available. "seekpath" modify the structure to put it into a conventional one.
  "legacy" is old method known to be wrong for ... but it doesn't modify the cell.  

* **protocol**, python str.
  Explained above


Therefore it is sufficient to have something like::

        bands_inp_gen = SiestaBandsInputsGenerator()
        builder = bands_inp_gen.get_builder(
                structure=structure, calc_engines=calc_engines, protocol=protocol, 
                path_generator=path_generator)


To obtain a builder that can be directly submitted (`submit(builder)`) and it performs the task of calculating
the electronic band structure of a system.
For the outputs, see specifications of the **SiestaBandsWorkchain**. Any input can be modified manually before
submission as the builder is not a stored object.
An example of the use is in aiida_siesta/examples/workflows/bands_with_protocol.py.


#Building nested workchains
#--------------------------
#
#One more type of example is presented in the folder aiida_siesta/examples/workflows/, the one with prefix
#`workchain_`. They show two simple example of nested workchains. The **SiestaBandsWorkchain** is launched
#inside another workchain, some post process of outputs is performed and new outputs returned.
