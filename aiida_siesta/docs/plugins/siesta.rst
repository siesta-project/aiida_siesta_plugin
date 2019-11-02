Standard Siesta plugin
++++++++++++++++++++++

Description
-----------

A plugin for Siesta's basic functionality. 


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1 series, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).

Inputs
------

* **code**, class :py:class:`Code <aiida.orm.Code>`

A code object linked to a Siesta executable. 
  
* **structure**, class :py:class:`StructureData <aiida.orm.StructureData>`

A structure. Siesta employs "species labels" to implement special
conditions (such as basis set characteristics) for specific atoms
(e.g., surface atoms might have a richer basis set). This is
implemented through the 'name' attribute of the Site objects. For example::

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


* **parameters**, class :py:class:`Dict <aiida.orm.Dict>`

A dictionary with scalar fdf variables and blocks, which are the
basic elements of any Siesta input file. A given Siesta fdf file
can be cast almost directly into this dictionary form, except that
some items (for structure data) are blocked. Any units are
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
    }

Note that Siesta fdf keywords allow '.', '-', (or nothing) as internal
separators. AiiDA does not allow the use of '.' in nodes to be
inserted in the database, so it should not be used in the input script
(or removed before assigning the dictionary to the Dict
instance). For legibility, a single dash ('-') is suggested, as in the
examples above.

* **pseudos**, input namespace of class :py:class:`PsfData  <aiida_siesta.data.psf.PsfData>`
  OR class :py:class:`PsmlData  <aiida_siesta.data.psml.PsmlData>`

The PsfData and PsmlData classes have been implemented along the lines of the Upf class for QE.

One pseudopotential file per atomic element. Several species (in the
Siesta sense, which allows the same element to be treated differently
according to its environment) can share the same pseudopotential. For the example
above::

  pseudo_file_to_species_map = [ ("C.psf", ['C', 'Cred']),
                            ("H.psf", 'H')
			    ]


Alternatively, a pseudo for every atomic species can be set with the
**use_pseudos_from_family**  method, if a family of pseudopotentials
has been installed. (But the family approach does not yet support
multiple species sharing the same pseudopotential.)

.. note:: The verdi command-line interface now supports entry points
   defined by external packages. We have implemented  `verdi data
   psf` and `verdi data psml` suites of commands: `uploadfamily`, `exportfamily`, and
   `listfamilies`.

It can be argued that a single `SiestaPseudo` class, with psf and psml
subclasses, might have been implemented. But the `PsmlData` class aims
to transcend Siesta and to be used by other codes.

* **basis**, class :py:class:`Dict  <aiida.orm.Dict>`
  
A dictionary specifically intended for basis set information. It
follows the same structure as the **parameters** element, including
the allowed use of fdf-block items. This raw interface allows a
direct translation of the myriad basis-set options supported by the
Siesta program. In future we might have a more structured input for
basis-set information.

* **kpoints**, class :py:class:`KpointsData <aiida.orm.KpointsData>`
  
Reciprocal space points for the full sampling of the BZ during the
self-consistent-field iteration. It must be given in mesh form. There is no support
yet for Siesta's kgrid-cutoff keyword.

If this node is not present, only the Gamma point is used for sampling.

* **bandskpoints**, class :py:class:`KpointsData <aiida.orm.KpointsData>`
  
Reciprocal space points for the calculation of bands.  They can be
given as a simple list of k-points, as segments with start and end
point and number of points, or as a complete automatic path, using the
functionality of modern versions of the class.

If this node is not present, no band structure is computed.

* **settings**, class  :py:class:`Dict <aiida.orm.Dict>`
      
An optional dictionary that activates non-default operations. For a list of possible
values to pass, see the section on :ref:`advanced features <siesta-advanced-features>`.

Outputs
-------

There are several output nodes that can be created by the plugin,
according to the calculation details.  All output nodes can be
accessed with the ``calculation.out`` method.

The output parser takes advantage of the structured output available
in Siesta as a Chemical Markup Language (CML) file. The CML-writer
functionality should be compiled in and active in the run!

* **output_parameters** :py:class:`Dict <aiida.orm.Dict>` 

A dictionary with metadata, scalar result values, a warnings
list, and possibly a timing section.
Units are specified by means of an extra item with '_units'
appended to the key::

    {
      "siesta:Version": "siesta-4.0.2",
      "E_fermi": -3.24,
      "E_fermi_units": "eV",
      "FreeE": -6656.2343
      "FreeE_units": "eV",
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

The scalar quantities to include are specified in a global-variable
in the parser. Currently they are the Kohn-Sham, Free, Band, and Fermi
energies, and the total spin. These are converted to 'float'.
As this dictionary is sorted, keys for program values and metadata are
intermixed.

The timing information (if present), includes the global walltime in
seconds, and a decomposition by sections of the code. Most relevant
are typically the `compute_DM` and `setup_H` sections.

The 'warnings' list contains program messages, labeled as INFO,
WARNING, or FATAL, read directly from a MESSAGES file produced by
Siesta, which include items from the execution of the program and
also a possible 'out of time' condition. This is implemented by
passing to the program the wallclock time specified in the script,
and checking at each scf step for the walltime consumed. This
'warnings' list can be examined by the parser itself to raise an
exception in the FATAL case.

* **forces_and_stress** :py:class:`ArrayData <aiida.orm.ArrayData>`

Contains the final forces (eV/Angstrom) and stresses (GPa) in array form.
  

* **output_structure** :py:class:`StructureData <aiida.orm.StructureData>`
  
Present only if the calculation is moving the ions.  Cell and ionic
positions refer to the last configuration.

* **bands**, :py:class:`BandsData  <aiida.orm.BandsData>`
  
Present only if a band calculation is requested (signaled by the
presence of a **bandskpoints** input node of class KpointsData)
Contains an array with the list of electronic energies for every
kpoint. For spin-polarized calculations, there is an extra dimension
for spin.
  
No trajectories have been implemented yet.

Errors
------

Errors during the parsing stage are reported in the log of the calculation (accessible 
with the ``verdi process report`` command). 
Moreover, they are stored in the `output_parameters` node under the key ``warnings``.

Restarts
--------

A restarting capability is implemented through the optional input
**parent_calc_folder**, :py:class:`RemoteData  <aiida.orm.RemoteData>`

which represents the remote scratch folder for a previous calculation.

The density-matrix file is copied from the old calculation scratch
folder to the new calculation's one.

This approach enables continuation of runs which have failed due to
lack of time or insufficient convergence in the allotted number of
steps.

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
