Base workflow
+++++++++++++

Description
-----------

The SIESTA program is able to perform, in a single run, the
computation of the electronic structure, the optional relaxation of
the input structure, and a final analysis step in which a variety of
magnitudes can be computed: band structures, projected densities of
states, etc. The operations to be carried out are specified in a very
flexible input format.  Accordingly, the **SiestaBaseWorkChain**
has been designed to be able to run the most general SIESTA
calculation, with support for most of the available options (limited
only by corresponding support in the parser plugin). The option specifications
of the **SiestaBaseWorkChain** follow the conventions already presented in the
:ref:`Siesta plugin <siesta-plugin-inputs>`. Therefore, for instance, the addition of
the input keyword **bandskpoints** triggers the calculation of the band structure
of a system, while it is sufficient to add the SIESTA MD keywords to the
**parameters** input in order to perforem the relaxation of a structure.
In contarst to the **SiestaCalculation** plugin, however, the 
workchain is able to automatically restart a calculation in case of failure (lack of
electronic-structure or geometry relaxation convergence, termination due to
walltime restrictions, etc).
Therefore, the **SiestaBaseWorkChain** is the suggested tool to run Siesta calculations
in the AiiDA framework. In fact, it retains the same level of flexibility of the most
general Siesta calculation, but it adds robusness thanks to its ability
to automatically respond to erros.
Examples on the use of the **SiestaBaseWorkChain** are presented in the folder
`/aiida_siesta/examples/workflows`.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the      
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.



.. _siesta-base-wc-inputs:

Inputs
------

Most inputs of the WorkChain are mirroring the siesta plugin inputs. Therefore, more
detailed information on them can be found :ref:`here <siesta-plugin-inputs>`.
The only difference is regarding the way the computational resources are passed.
The siesta plugin make use of ``metadada.options`` for this task, here, instead, we have
a dedicated input node. This node is the first point in the following list, describing
all the inputs of the WorkChain.

* **options**, class :py:class:`Dict <aiida.orm.Dict>`, *Mandatory*

  Execution options. In this dictionary the computational resources and
  scheduler specifications (queue, account, etc ..) must be specified.
  An example is::
        options = Dict(
             dict={
                'max_wallclock_seconds': 360,
                'withmpi': True,
                'account': 'tcphy113c',
                'queue_name': 'DevQ',
                'resources': {'num_machines': 1,'num_mpiprocs_per_machine': 2},
                }
             )

  The `resources` and `max_wallclock_seconds` are required by AiiDA, the rest of the options
  depend on the scheduler of the machine one is submitting to.

.. |br| raw:: html

    <br />

* **code**,  class :py:class:`Code  <aiida.orm.Code>`, *Mandatory*

  A database object representing a Siesta executable. See the plugin documentation for more details.

.. |br| raw:: html

    <br />

* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`, *Mandatory*

  A structure. See the plugin documentation for more details.

.. |br| raw:: html

    <br />

* **parameters**, class :py:class:`Dict <aiida.orm.Dict>`,  *Mandatory*

  A dictionary with scalar fdf variables and blocks, which are the
  basic elements of any Siesta input file. A given Siesta fdf file
  can be cast almost directly into this dictionary form, except that
  some items (e.g. for structure data) are blocked. Any units are
  specified for now as part of the value string. Blocks are entered
  by using an appropriate key and Python's multiline string
  constructor. For example::
  
      {
        "mesh-cutoff": "200 Ry",
        "dm-tolerance": "0.0001",
        "%block example-block":
  	  """
  	  first line
  	  second line             
  	  %endblock example-block""",
        ...
      }
  
  Note that Siesta fdf keywords allow '.', '-', or nothing as
  internal separators. AiiDA does not allow the use of '.' in
  nodes to be inserted in the database, so it should not be used
  in the input script (or removed before assigning the dictionary to
  the Dict instance). For legibility, a single dash ('-') is suggested, as in the
  examples above. See the plugin documentation for more details on the blocked
  items.

.. |br| raw:: html

    <br />

* **pseudos**, input namespace of class :py:class:`PsfData <aiida_siesta.data.psf.PsfData>`
  OR class :py:class:`PsmlData <aiida_siesta.data.psml.PsmlData>`, *Optional*

  A dictionary of `PsfData  <aiida_siesta.data.psf.PsfData>` or
  `PsmlData  <aiida_siesta.data.psml.PsmlData>` objects representing the pseudopotentials for
  the calculation. See the plugin documentation for more details.
  In contrast to the case of the siesta plugin, the **pseudos** input
  is not mandatory. The **SiestaBaseWorkChain** supports, in fact, the direct use of
  **pseudo_family** (see below). If **pseudos** is not in input, a **pseudo_family** 
  specification must be used.

.. |br| raw:: html

    <br />

* **pseudo_family**, class :py:class:`Str <aiida.orm.Str>`, *Optional*

  String representing the name of a pseudopotential family stored in the database.
  Pseudofamilies can be uploaded in the database via the ``verdi data psf uploadfamily``
  or ``verdi data psml uploadfamily`` CLI interface.

.. |br| raw:: html

    <br />

* **basis**, class :py:class:`Dict  <aiida.orm.Dict>`, *Optional*
  
  A dictionary specifically intended for basis set information. It
  follows the same structure as the **parameters** element, including
  the allowed use of fdf-block items. This raw interface allows a
  direct translation of the myriad basis-set options supported by the
  Siesta program. If not specified, a calculation with only the gamma 
  point is performed. See the plugin documentation for more details.

.. |br| raw:: html

    <br />

* **kpoints**, class :py:class:`KpointsData <aiida.orm.KpointsData>`, *Optional*
  
  Reciprocal space points for the full sampling of the BZ during the
  self-consistent-field iteration. It must be given in mesh form. There is no support
  yet for Siesta's kgrid-cutoff keyword. See the plugin documentation for more details.
  If this node is not present, only the Gamma point is used for sampling.

.. |br| raw:: html

    <br />

* **bandskpoints**, class :py:class:`KpointsData  <aiida.orm.KpointsData>`, *Optional*
  
  Reciprocal space points for the calculation of bands.  They can be
  given as a simple list of k-points, as segments with start and end
  point and number of points, or as a complete automatic path, using the
  functionality of modern versions of the class. See the plugin documentation 
  for more details.
  If this node is not present, no band structure is computed.

.. |br| raw:: html

    <br />

* **settings**, class :py:class:`Dict <aiida.orm.Dict>`, *Optional*
      
  An optional dictionary that activates non-default operations. For a list of possible
  values to pass, see the section on :ref:`advanced features <siesta-advanced-features>`.

.. |br| raw:: html

    <br />

* **clean_workdir**, class :py:class:`Bool <aiida.orm.Bool>`, *Optional*

  If true, work directories of all called calculations will be cleaned
  out. Default is false.

.. |br| raw:: html

    <br />

* **max_iterations**, class :py:class:`Int <aiida.orm.Int>`, *Optional*

  The maximum number of iterations allowed in the restart cycle for
  calculations. The **SiestaBaseWorkChain** tries to deal with some 
  common siesta errors (see `here <basewc-error>`) and restart the calculation with appropriate
  modifications. The integer **max_iterations** is the maximum number
  of times the restart is performed no matter what error is recorded.
  The input is optional, if not specified, the default `Int(5)` is used.

.. |br| raw:: html

    <br />

* **parent_calc_folder**, class  :py:class:`RemoteData <aiida.orm.RemoteData>` , *Optional*

  Optional port used to activate the restart features, as explained in the plugin documentation.


Relaxation and bands
--------------------
As already mentioned in the introduction, in addition to simple scf calculations, the **SiestaBaseWorkChain** 
can be used to perform the relaxation of a structure and the electronic bands calculations.
For the electronic bands, however, we suggest the use of the **BandgapWorkChain** distributed in this package, because
it adds the feature to automatically calculate the band gap.
Concerning the relaxation of a structure, the **SiestaBaseWorkChain** simply exploits the internal relaxation
implemented in Siesta in order to complete the task. The full set of a Siesta relaxation options can be
accessed just adding the corresponding keyword and value in the **parameters** input dictionary. The only additional
feature that the **SiestaBaseWorkChain** adds is that it requires to reach the target forces and stress
to consider completed the task. If this does not happen in a single Siesta run, the workchain restarts
automatically the relaxation. The maximum number of restarts is specified with the keyword **max_iterations**,
as explained in the previous subsection.


Submitting the WorkChain
------------------------

WorkChains are submitted in AiiDA exacly like any other calculation. Therefore::

        from aiida_siesta.workflows.base import SiestaBaseWorkChain
        from aiida.engine import
        builder = SiestaBaseWorkChain.get_builder()
        builder.options = options
        ... All the inputs here ...
        submit(builder) #or run

There is no need to set the computational resources with the metadata as they are already
defined in the input **options**, however ``builder.metadata.label`` and ``builder.metadata.description``
could be used to label and describe the WorkChain.
Again, the use of the ``builder`` is not mandatory, the inputs can be passed as arguments of
``sumbit``/``run`` as explained in the siesta plugin section.

Outputs
-------

The outputs of the **SiestaBaseWorkChain** mirror exactly the one of the siesta plugin.
Therefore all the information can be obtained in the corresponding section.
We list here the outputs.

* **output_parameters** :py:class:`Dict <aiida.orm.Dict>` 

  A dictionary with metadata and scalar result values from the last
  calculation executed.

.. |br| raw:: html

    <br />

* **output_structure** :py:class:`StructureData <aiida.orm.StructureData>`
  
  Present only if the workchain is modifying the geometry of the system.

.. |br| raw:: html

    <br />

* **bands**, :py:class:`BandsData <aiida.orm.BandsData>`
  
  Present only if a band calculation is requested (signaled by the
  presence of a **bandskpoints** input node of class KpointsData)
  Contains an array with the list of electronic energies for every
  kpoint. For spin-polarized calculations, there is an extra dimension
  for spin.

.. |br| raw:: html

    <br />

* **forces_and_stress** :py:class:`ArrayData <aiida.orm.ArrayData>`

  Contains the final forces (`eV/Angstrom`) and stresses (`GPa`) in array form.

.. |br| raw:: html

    <br />

* **remote_folder**, :py:class:`RemoteData <aiida.orm.RemoteData>`

  The working remote folder for the last calculation executed. As the **SiestaBaseWorkChain**
  automatically restarts the calculation in case of common failures, the very last
  siesta calculation is considered the interesting one for a further manual restart.
  Therefore its folder is returned in this node.


.. _basewc-error:

Error handling
--------------

We list here the errors that are handled by the **SiestaBaseWorkChain** and the
corresponding action taken. The error are actually detected by the siesta parser,
in the WorkChain, the handling is performed.

* **SCF_NOT_CONV**

  When the convergence of the self-consistent cycle is not reached in ``max-scf-iterations`` or
  in the allocated ``max_walltime``, siesta raises the **SCF_NOT_CONV** error.
  The **SiestaBaseWorkChain** is able to detect this error and restart the calculation with no
  modifications on the input parameters.

.. |br| raw:: html

    <br />

* **GEOM_NOT_CONV**

  When the convergence of the geometry (during a relaxation) is not reached
  in the allocated ``max_walltime``, siesta raises the **GEOM_NOT_CONV** error.
  The **SiestaBaseWorkChain** is able to detect this error and restart the calculation with no
  modifications on the input parameters.

.. |br| raw:: html

    <br />

* **SPLIT_NORM**

  The **SiestaBaseWorkChain** deals with problems connected to the basis set creation.
  If a "too small split-norm" error is detected, the WorkChains reacts in two ways.
  If a global split-norm was defined in input through ``pao-split-norm``, its value is reset to
  the minimum acceptable. If no global split-norm was defined the option ``pao-split-tail-norm = True``
  is set.

Two more errors are detected by the WorkChain, but not handled at the moment,
only a specific error code is returned as output without attempting a restart.

* **BASIS_POLARIZ**

  If an error on the polarization of one orbital is detected, the error code 403 is returned.
  The solution to this problem is to set the "non-perturbative" polarization scheme for the
  element that presents an error, however this possibility is available only in recent
  versions of AiiDA, making inconvenient to treat automatically the resolution of this error.

.. |br| raw:: html

    <br />

* **ERROR_BANDS**

  If a calculation of the electronic bands is requested, but
  an error in the parsing of the bands file is detected, the error code 404 is returned.
  In this case, the WorkChain will anyway return all the other outputs because the checks
  on the bands file are always performed at the very end of the calculation.

The **SiestaBaseWorkChain** also inherits the error codes of the **BaseRestartWorkChain**
of the aiida-core distribution. For instance,
if an unexpected error is raised twice, the workchain finishes with exit code 402, if the
maximum number of iterations is reached, error 401 is returned. More in the section
`BaseRestartWorkChain`_ of the aiida-core package.

Protocol system
---------------

The protocol system is available for this WorkChain. The ``SiestaBaseWorkChain.inputs_generator()``
makes available all the methods explained in the :ref:`protocols documentation <how-to>`. For example::

        from aiida_siesta.workflows.base import SiestaBaseWorkChain
        inp_gen = SiestaBaseWorkChain.inputs_generator()
        builder = inp_gen.get_filled_builder(structure, calc_engines, protocol)
        #here user can modify builder befor submission.
        submit(builder)

is sufficient to submit a **SiestaBaseWorkChain** on ``structure`` following the specifications of
``protocols`` and computational resources collected in ``calc_engines``.
The structure of ``calc_engines`` is the same as for the **SiestaCalculation** input generator
(again see :ref:`protocols documentation <how-to>`).


.. _BaseRestartWorkChain: https://aiida.readthedocs.io/projects/aiida-core/en/latest/reference/apidoc/aiida.engine.processes.html?highlight=baserestart#aiida.engine.processes.BaseRestartWorkChain
