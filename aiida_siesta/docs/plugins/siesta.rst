Siesta calculations
+++++++++++++++++++

Description
-----------

A plugin for Siesta main code. It allows to prepare, submit and retrieve the results of a standard siesta calculation,
including support for the parsing of the electronic bands and the output geometry of a relaxation. 
It is implemented in the class **SiestaCalculation**.

Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform (https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the 
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.

.. _siesta-plugin-inputs:

Inputs
------
Some examples are referenced in the following list. They are located in the folder
`aiida_siesta/examples/plugins/siesta`.

* **code**, class :py:class:`Code <aiida.orm.Code>`, *Mandatory*

  A code object linked to a Siesta executable. 
  If you setup the code ``Siesta4.0.1`` on machine ``kelvin`` following the `aiida guidelines`_,
  then the code is selected in this way::

        codename = 'Siesta4.0.1@kelvin'
        from aiida.orm import Code
        code = Code.get_from_string(codename)

 
* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`, *Mandatory*

  A structure. Siesta employs "species labels" to implement special
  conditions (such as basis set characteristics) for specific atoms
  (e.g., surface atoms might have a richer basis set). This is
  implemented through the ``name`` attribute of the Site objects. For example::
      from aiida.orm import StructureData

      alat = 15. # angstrom
      cell = [[alat, 0., 0.,],
        [0., alat, 0.,],
        [0., 0., alat,],
       ]

       # Benzene molecule with a special carbon atom
       s = StructureData(cell=cell)
       s.append_atom(position=(0.000,0.000,0.468),symbols=['H'])
       s.append_atom(position=(0.000,0.000,1.620),symbols=['C'])
       s.append_atom(position=(0.000,-2.233,1.754),symbols=['H'])
       s.append_atom(position=(0.000,2.233,1.754),symbols=['H'])
       s.append_atom(position=(0.000,-1.225,2.327),symbols='C',name="Cred")
       s.append_atom(position=(0.000,1.225,2.327),symbols=['C'])
       s.append_atom(position=(0.000,-1.225,3.737),symbols=['C'])
       s.append_atom(position=(0.000,1.225,3.737),symbols=['C'])
       s.append_atom(position=(0.000,-2.233,4.311),symbols=['H'])
       s.append_atom(position=(0.000,2.233,4.311),symbols=['H'])
       s.append_atom(position=(0.000,0.000,4.442),symbols=['C'])
       s.append_atom(position=(0.000,0.000,5.604),symbols=['H'])

  The class :py:class:`StructureData <aiida.orm.StructureData>` uses Angstrom
  as internal units, the cell and atom positions must be specified in Angstrom.

  The :py:class:`StructureData <aiida.orm.StructureData>` can also import 
  ase structures or pymatgen structures. These two tools can be used to load
  structure from files. See example `example_cif_bands.py`.

.. |br| raw:: html

    <br />

* **parameters**, class :py:class:`Dict <aiida.orm.Dict>`, *Mandatory*

  A dictionary with scalar fdf variables and blocks, which are the
  basic elements of any Siesta input file. A given Siesta fdf file
  can be cast almost directly into this dictionary form, except that
  some items are blocked. The blocked keywords include the system information
  (``system-label``, ``system-name``) and all the structure information as they
  will be automatically set by Aiida. Moreover, the keyword ``dm-use-save-dm`` is
  not allowed (the restart options are explained :ref:`here <siesta-restart>`)
  together with the keyword ``geometry-must-converge`` (set to True by default for each
  calculation with variable geometry). Finally,  all the ``pao`` options must be avoided here, 
  because they belong to the **basis** input 
  (next to next in this list). Any units are
  specified for now as part of the value string. Blocks are entered
  by using an appropriate key and Python's multiline string
  constructor. For example::

    from aiida.orm import Dict    

    parameters = Dict(dict={
      "mesh-cutoff": "200 Ry",
      "dm-tolerance": "0.0001",
      "%block example-block":
        """
        first line
        second line
        %endblock example-block""",
    })

  Note that Siesta fdf keywords allow '.', '-', (or nothing) as internal
  separators. AiiDA does not allow the use of '.' in nodes to be
  inserted in the database, so it should not be used in the input script
  (or removed before assigning the dictionary to the Dict
  instance). For legibility, a single dash ('-') is suggested, as in the
  examples above. Moreover, because the parameters are passed through a python 
  dictionary, if, by mistake, the user passes the same keyword two (or more) 
  times, only the last specification will be considered. For instance::

     parameters = Dict(dict={
       "mesh-cutoff": "200 Ry",
       "mesh-cutoff": "300 Ry",
       })

  will set a ``mesh-cutoff`` of `300 Ry`. This is the opposite respect to what is done 
  in the Siesta code, where the first assignment is the selected one. Please note that 
  this applies also to keywords that correspond to the same fdf variable. For instance::

     parameters = Dict(dict={
       "mesh-cutoff": "200 Ry",
       "Mesh-Cut-off": "300 Ry",
       })

  will run a calculation with ``mesh-cutoff`` equal to `300 Ry`, whithout raising any
  error.


.. |br| raw:: html

    <br />

* **pseudos**, input namespace of class :py:class:`PsfData  <aiida_siesta.data.psf.PsfData>`
  OR class :py:class:`PsmlData  <aiida_siesta.data.psml.PsmlData>`, *Mandatory*

  The `PsfData  <aiida_siesta.data.psf.PsfData>` and `PsmlData  <aiida_siesta.data.psml.PsmlData>`
  classes have been implemented along the lines of the Upf class of aiida-core.

  One pseudopotential file per atomic element is required. Several species (in the
  Siesta sense, which allows the same element to be treated differently
  according to its environment) can share the same pseudopotential. For the example
  above::

    import os
    from aiida_siesta.data.psf import PsfData

    pseudo_file_to_species_map = [ ("C.psf", ['C', 'Cred']),("H.psf", ['H'])]
    pseudos_dict = {}
    for fname, kinds, in pseudo_file_to_species_map:
          absname = os.path.realpath(os.path.join("path/to/file",fname))
          pseudo, created = PsfData.get_or_create(absname, use_first=True)
          for j in kinds:
                  pseudos_dict[j]=pseudo

  Alternatively, a pseudo for every atomic species can be set with the
  ``use_pseudos_from_family``  method, if a family of pseudopotentials
  has been installed. For an example see `example_psf_family.py`

  .. note:: The verdi command-line interface now supports entry points
     defined by external packages. We have implemented  `verdi data
     psf` and `verdi data psml` suites of commands: `uploadfamily`, `exportfamily`, and
     `listfamilies`.

  It can be argued that a single "SiestaPseudo" class, with psf and psml
  subclasses, might have been implemented. But the `PsmlData  <aiida_siesta.data.psml.PsmlData>`
  class aims to transcend Siesta and to be used by other codes.

.. |br| raw:: html

    <br />

* **basis**, class :py:class:`Dict  <aiida.orm.Dict>`, *Optional*

  A dictionary specifically intended for basis set information. It
  follows the same structure as the **parameters** element, including
  the allowed use of fdf-block items. This raw interface allows a
  direct translation of the myriad basis-set options supported by the
  Siesta program. In future we might have a more structured input for
  basis-set information.
  An example::

        from aiida.orm import Dict

        basis_dict = {
        'pao-basistype':'split',
        'pao-splitnorm': 0.150,
        'pao-energyshift': '0.020 Ry',
        '%block pao-basis-sizes':
        """
        C    SZP
        Cred SZ
        H    SZP
        %endblock pao-basis-sizes""",
        }

        basis = Dict(dict=basis_dict)

  In case no basis is set, the Siesta calculation will not include
  any basis specification and it will run with the default basis: DZP 
  plus (many) other defaults.

.. |br| raw:: html

    <br />

* **kpoints**, class :py:class:`KpointsData <aiida.orm.KpointsData>`, *Optional*

  Reciprocal space points for the full sampling of the BZ during the
  self-consistent-field iteration. It must be given in mesh form. There is no support
  yet for Siesta's "kgrid-cutoff" keyword::
          from aiida.orm import KpointsData
          kpoints=KpointsData()
          kp_mesh = 5
          mesh_displ = 0.5 #optional
          kpoints.set_kpoints_mesh([kp_mesh,kp_mesh,kp_mesh],[mesh_displ,mesh_displ,mesh_displ])

  The class `KpointsData <aiida.orm.KpointsData>` also implements the methods 
  ``set_cell_from_structure`` and ``set_kpoints_mesh_from_density``
  that allow to obtain a uniform mesh automatically.
  
  If this node is not present, only the Gamma point is used for sampling.

.. |br| raw:: html

    <br />

* **bandskpoints**, class :py:class:`KpointsData <aiida.orm.KpointsData>`, *Optional*

  Reciprocal space points for the calculation of bands.
  The **full list of kpoints must be passed** to ``bandskpoints``
  and they must be in **units of the reciprocal lattice vectors**.
  There is no obligation to set the cell in ``bandskpoints``, however this might be useful
  in order to exploit the functionality
  of the class :py:class:`KpointsData <aiida.orm.KpointsData>`.
  If set, the cell must be the same of the input **structure**.
  Some examples on how to pass the kpoints are the following.

  One can manually listing a set of isolated kpoints::
          from aiida.orm import KpointsData
          bandskpoints=KpointsData()
          kpp = [(0.1,  0.1, 0.1), (0.5,  0.5, 0.5), (0., 0., 0.)]
          bandskpoints.set_kpoints(kpp)
  In this case the Siesta input will use the "BandPoints" block.
  
  Alternatively (recommended) the high-symmetry path associated to the
  structure under investigation can be
  automatically generated through the aiida tool ``get_explicit_kpoints_path``.
  Here how to use it::
          from aiida.orm import KpointsData
          bandskpoints=KpointsData()
          from aiida.tools import get_explicit_kpoints_path
          symmpath_parameters = Dict(dict={'reference_distance': 0.02})
          kpresult = get_explicit_kpoints_path(s, **symmpath_parameters.get_dict())
          bandskpoints = kpresult['explicit_kpoints']
  Where 's' in the input structure and ``reference_distance`` is
  the distance between two subsequent kpoints.
  In this case the block "BandLines" is set in the Siesta
  calculation.

  .. warning:: "SeeK-path"
     might modify the structure to follow particular conventions
     and the generated kpoints might only 
     apply on the internally generated 'primitive_structure' and not 
     on the input structure that was provided. The correct
     way to use this tool is to use the generated 'primitive_structure' also for the
     Siesta calculation::
          structure = kpresult['primitive_structure']

  .. warning:: As we use the initial structure cell in order to obtain
     the kpoints path, it is very risky to apply this method when also a relaxation
     of the cell is performed!
     The cell might relax in a different symmetry resulting in a wrong
     path for the bands.
     Consider to use the `BandGapWorkChain` if a relaxation is needed
     before computing the bands.

  .. note:: The ``get_explicit_kpoints_path`` make use of "SeeK-path".
     Please cite the `HPKOT paper`_ if you use this tool. "SeeK-path"
     is a external utility, not a requirement for aiida-core, therefore
     it is not available by default. It can be easily installed using
     ``pip install seekpath``. "SeeK-path" allows to
     determine canonical unit cells and k-point information in an easy
     way. For more general information, refer to the `SeeK-path documentation`_.


  The final option covers the situation
  when one needs to calculate the bands on a specific path
  (and maybe needs to maintain a specific convention for the
  structure). The full list of kpoints must be passed and, very
  importantly, labels must be set for the high symmetry points!
  This is essential for the correct set up of the "BandLines" in Siesta.
  External tolls can be used to create equidistant points, whithin aiida
  the following (very involved) option is available::
        from aiida.orm import KpointsData
        bandskpoints=KpointsData()
        from aiida.tools.data.array.kpoints.legacy import get_explicit_kpoints_path as legacy_path
        kpp = [('A',  (0.500,  0.250, 0.750), 'B', (0.500,  0.500, 0.500), 40),
        ('B', (0.500,  0.500, 0.500), 'C', (0., 0., 0.), 40)]
        tmp=legacy_path(kpp)
        bandskpoints.set_kpoints(tmp[3])
        bandskpoints.labels=tmp[4]
  The legacy ``get_explicit_kpoints_path`` shares only the name with the function in
  ``aiida.tools``, but it is very different in scope.

  The full list of cases can be explored looking at the example example_bands.py

  .. warning:: The implementation relies on the correct description of
     the labels in the class :py:class:`KpointsData <aiida.orm.KpointsData>`.
     Refrain from improper use of ``bandskpoints.labels`` and follow the 
     the instructions described above. An incorrect use of the labels
     might result in an incorrect parsing of the bands.

  If the keyword node **bandskpoints** is not present, no band structure is computed.

.. |br| raw:: html

    <br />

* **settings**, class  :py:class:`Dict <aiida.orm.Dict>` , *Optional*      

  An optional dictionary that activates non-default operations. For a list of possible
  values to pass, see the section on :ref:`advanced features <siesta-advanced-features>`.

.. |br| raw:: html

    <br />

* **parent_calc_folder**, class  :py:class:`RemoteData <aiida.orm.RemoteData>` , *Optional*

  Optional port used to activate the :ref:`restart features <siesta-restart>`.


Submitting the calculation
--------------------------

Once all the inputs above are set, the subsequent step consists in passing them to the 
calculation class and run/submit it.

First, the Siesta calculation class is loaded::

        from aiida_siesta.calculations.siesta import SiestaCalculation
        builder = SiestaCalculation.get_builder()

The inputs (defined as in the previous section) are passed to the builder::

        builder.code = code
        builder.structure = structure
        builder.parameters = parameters
        builder.pseudos = pseudos_dict
        builder.basis = basis
        builder.kpoints = kpoints
        builder.bandskpoints = bandskpoints

Finally the resources for the calculation must be set, for instance::

        builder.metadata.options.resources = {'num_machines': 1}
        builder.metadata.options.max_wallclock_seconds = 1800

Optionally, label and description::

        builder.metadata.label = 'My generic title'
        builder.metadata.description 'My more detailed description'

To run the calculation in an interactive way::

        from aiida.engine import run
        results = run(builder)
Here the results variable will contain a dictionary 
containing all the nodes that were produced as output.

Another option is to submit it to the daemon::

        from aiida.engine import submit
        calc = submit(builder)
In this case, calc is the calculation node and not the results dictionary.

.. note:: In order to inspect the inputs created by AiiDA without 
   actually running the calculation, we can perform a dry run of the submission process::
        builder.metadata.dry_run = True
        builder.metadata.store_provenance = False
   This will create the input files, that are available for inspection.

.. note:: The use of the builder makes the process more intuitive, but it
   is not mandatory. The inputs can be provided as keywords argument when you 
   launch the calculation, passing the calculation class as the first argument::
        run(SiestaCalculation, structure=s, pseudos=pseudos, kpoints = kpoints, ...)
   same syntax for the command ``submit``.

A large set of examples covering some standard cases are in the folder 
`aiida_siesta/examples/plugins/siesta`. They can be run with::
        runaiida example_name.py {--send, --dont-send} code@computer

The parameter ``--dont-send`` will activate the "dry run" option. In that case a test
folder (`submit_test`) will be created, containing all the files that aiida
generates automatically. The parameter ``--send`` will submit the example to
the daemon. One of the two options needs to be present to run the script. 
The second argument contains the name of the code (``code@computer``) to use
in the calculation. It must be a previously set up code, corresponding to
a siesta executable.


Outputs
-------

There are several output nodes that can be created by the plugin,
according to the calculation details.  All output nodes can be
accessed with the ``calculation.outputs`` method.


* **output_parameters** :py:class:`Dict <aiida.orm.Dict>` 

  A dictionary with metadata, scalar result values, a warnings
  list, and possibly a timing section.
  Units are specified by means of an extra item with '_units'
  appended to the key::
    {
      "siesta:Version": "siesta-4.0.2",
      "E_Fermi": -3.24,
      "E_Fermi_units": "eV",
      "FreeE": -6656.2343,
      "FreeE_units": "eV",
      "E_KS": -6656.2343,
      "E_KS_units": 'eV',
      "global_time": 55.213,
      "timing_decomposition": {
        "compute_DM": 33.208, 
        "nlefsm-1": 0.582, 
        "nlefsm-2": 0.045, 
        "post-SCF": 2.556, 
        "setup_H": 16.531, 
        "setup_H0": 2.351, 
        "siesta": 55.213, 
        "state_init": 0.171
      }, 
      "warnings": [ "INFO: Job Completed"]
    }

  The scalar quantities included are, currently, the Kohn-Sham
  (``E_KS``), Free (``FreeE``), Band (``Ebs``), and Fermi (``E_Fermi``)
  energies, and the total spin (``stot``). These are converted to :py:class:`float <float>`.
  The other quantities are or type :py:class:`str <str>`.

  The timing information (if present), includes the global walltime in
  seconds, and a decomposition by sections of the code. Most relevant
  are typically the ``compute_DM`` and ``setup_H`` sections.

  The ``warnings`` list contains program messages, labeled as "INFO",
  "WARNING", or "FATAL", read directly from a  `MESSAGES` file produced by
  Siesta, which include items from the execution of the program and
  also a possible 'out of time' condition. This is implemented by
  passing to the program the wallclock time specified in the script,
  and checking at each scf step for the walltime consumed. This
  ``warnings`` list can be examined by the parser itself to raise an
  exception in the "FATAL" case.

.. |br| raw:: html

    <br />

* **forces_and_stress** :py:class:`ArrayData <aiida.orm.ArrayData>`

  Contains the final forces (`eV/Angstrom`) and stresses (`GPa`) in array form.
  To access their values::
        forces_and_stress.get_array("forces")
        forces_and_stress.get_array("stress")
  
.. |br| raw:: html

    <br />

* **output_structure** :py:class:`StructureData <aiida.orm.StructureData>`
  
  Present only if the calculation is moving the ions.  Cell and ionic
  positions refer to the last configuration.

.. |br| raw:: html

    <br />

* **bands**, :py:class:`BandsData  <aiida.orm.BandsData>`
  
  Present only if a band calculation is requested (signaled by the
  presence of a **bandskpoints** input node of class `KpointsData <aiida.orm.KpointsData>`).
  It contains an array with the list of electronic energies (in `eV`) for every
  kpoint. For spin-polarized calculations, there is an extra dimension
  for spin. In this class also the full list of kpoints is stored and they are
  in units of `1/Angstrom`. Therefore a direct comparison with the Siesta output 
  SystLabel.bands is possible only after the conversion of `Angstrom` to `Bohr`.
  The bands are not rescaled by the Fermi energy. Tools for the generation
  of files that can be easly plot are available through ``bands.export``.

.. |br| raw:: html

    <br />

* **remote_folder**, :py:class:`RemoteData <aiida.orm.RemoteData>`

  The working remote folder for the last calculation executed.


.. |br| raw:: html

    <br />

* **retrieved**, :py:class:`RemoteData <aiida.orm.RemoteData>`

  The local folder with the retrieved files.


No trajectories have been implemented yet.

Errors
------

Errors during the parsing stage are reported in the log of the calculation (accessible 
with the ``verdi process report`` command). 
Moreover, they are stored in the **output_parameters** node under the key ``warnings``.

.. _siesta-restart:

Restarts
--------

A restarting capability is implemented through the optional input
**parent_calc_folder**, :py:class:`RemoteData  <aiida.orm.RemoteData>`,
which represents the remote scratch folder (**remote_folder** output)
of a previous calculation.

The density-matrix file is copied from the old calculation scratch
folder to the new calculation's one.

This approach enables continuation of runs which have failed due to
lack of time or insufficient convergence in the allotted number of
steps.

An informative example is `example_restart.py` in the folder `aiida_siesta/examples/plugins/siesta`.

.. _siesta-advanced-features:

Additional advanced features
----------------------------

While the input link with name **parameters** is used for the main
Siesta options (as would be given in an fdf file), additional settings
can be specified in the **settings** input, also of type Dict.

Below we summarise some of the options that you can specify, and their effect.

The keys of the settings dictionary are internally converted to
uppercase by the plugin.

Adding command-line options
...........................

If you want to add command-line options to the executable (particularly 
relevant e.g. to tune the parallelization level), you can pass each option 
as a string in a list, as follows::

  settings_dict = {  
      'cmdline': ['-option1', '-option2'],
  }
  builder.settings = Dict(dict=settings_dict)

Note that very few user-level comand-line options (besides those
already inserted by AiiDA for MPI operation) are currently implemented.

Retrieving more files
.....................

If you know that your calculation is producing additional files that you want to
retrieve (and preserve in the AiiDA repository), you can add
those files as a list as follows::


  settings_dict = {  
    'additional_retrieve_list': ['aiida.EIG', 'aiida.ORB_INDX'],
  }
   builder.settings = Dict(dict=settings_dict)

See for example `example_ldos.py` in `aiida_siesta/examples/plugins/siesta`.
The files can then be accesed through the output **retrieved** and
its methods ``get_object`` and ``get_object_content``.

.. _SeeK-path documentation: https://seekpath.readthedocs.io/en/latest/
.. _aiida guidelines: https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html
.. _HPKOT paper: http://dx.doi.org/10.1016/j.commatsci.2016.10.015
