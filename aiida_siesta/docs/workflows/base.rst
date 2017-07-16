SIESTA Base workflow
++++++++++++++++++++++

Description
-----------

The SIESTA program is able to perform, in a single run, the
computation of the electronic structure, the optional relaxation of
the input structure, and a final analysis step in which a variety of
magnitudes can be computed: band structures, projected densities of
states, etc. The operations to be carried out are specified in a very
flexible input format.  Accordingly, the **SiestaBaseWorkchain**
has been designed to be able to run the most general SIESTA
calculation, with support for most of the available options (limited
only by corresponding support in the parser plugin). In addition, the
workchain is able to restart a calculation in case of failure (lack of
electronic-structure or relaxation convergence, termination due to
walltime restrictions, etc).


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1 series, which
can be found in the development platform
(http://launchpad.net/siesta/).

Inputs
------

* **code**, a code

* **structure**, class :py:class:`StructureData
  <aiida.orm.data.structure.StructureData>`

A structure. See the plugin documentation for more details.


* **parameters**, class :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`

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
	  "%block example-block": """
	  first line
	  second line             """,
    }

Note that Siesta fdf keywords allow '.', '-', or nothing as
internal separators. AiiDA does not allow the use of '.' in
nodes to be inserted in the database, so it should not be used
in the input script (or removed before assigning the dictionary to
the ParameterData instance).

* **pseudos**

(Optional)
A dictionary of PsfData objects representing the pseudopotentials for
the calculation. If it is not input, a **pseudo_family** specification
must be used (see below).

The PsfData class has been implemented along the lines of the Upf class for QE.


* **pseudo_family**, a Str

(Optional)


* **basis**, class :py:class:`ParameterData  <aiida.orm.data.parameter.ParameterData>`
  
A dictionary specifically intended for basis set information. It
follows the same structure as the **parameters** element, including
the allowed use of fdf-block items. This raw interface allows a
direct translation of the myriad basis-set options supported by the
Siesta program. In future we might have a more structured input for
basis-set information.

* **kpoints**, class :py:class:`KpointsData <aiida.orm.data.array.kpoints.KpointsData>`
  
Reciprocal space points for the full sampling of the BZ during the
self-consistent-field iteration. It must be given in mesh form. There is no support
yet for Siesta's kgrid-cutoff keyword.

If this node is not present, only the Gamma point is used for sampling.

* **bandskpoints**, class :py:class:`KpointsData
  <aiida.orm.data.array.kpoints.KpointsData>`
  
Reciprocal space points for the calculation of bands.  They can be
given as a simple list of k-points, as segments with start and end
point and number of points, or as a complete automatic path, using the
functionality of modern versions of the class.

If this node is not present, no band structure is computed.

* **settings**, class
  :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`
      
An optional dictionary that activates non-default operations. For a list of possible
values to pass, see the section on :ref:`advanced features <siesta-advanced-features>`.

* **options**, class
  :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`

Execution options

* **clean_workdir**, Bool

* **max_iterations**, Int

The maximum number of iterations allowed in the restart cycle for
calculations.


Outputs
-------

* **output_parameters** :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>` 
  (accessed by ``calculation.res``)

A dictionary with metadata and scalar result values from the last
calculation executed.

* **output_structure** :py:class:`StructureData
  <aiida.orm.data.structure.StructureData>`
  
Present only if the workchain is modifying the geometry of the system.

* **bands_array**, :py:class:`BandsData
  <aiida.orm.data.array.bands.BandsData>`
  
Present only if a band calculation is requested (signaled by the
presence of a **bandskpoints** input node of class KpointsData)
Contains the list of electronic energies for every kpoint. For
spin-polarized calculations, the 'bands' array has an extra dimension
for spin.

* **remote_folder**

The working remote folder for the last calculation executed.


