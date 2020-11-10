Bandgap workflow
++++++++++++++++

Description
-----------

The **BandgapWorkChain** is an extension of the **SietaBaseWorkChain** 
that introduces some logic to automatically obtain the bands and
applyes a simple post-process with the scope to return the metallic or
insulating nature of the material and, possibly, the band gap.

To calculate the gap, this workchain makes use of a tool distributed in aiida-core,
the method ``find_bandgap`` hosted in ``aiida.orm.nodes.data.array.bands``.

The optional automatic generation of the kpoints path for the bands
is done using `SeeK-path`_.

Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.


Inputs
------

All the **SiestaBaseWorkChain** inputs are as well inputs of the **BangapWorkChain**,
therefore the system and DFT specifications (structure, parameters, etc.) are
inputted in the WorkChain using the same syntax explained in the **SiestaBaseWorkChain**
:ref:`documentation <siesta-base-wc-inputs>`.
There is however the addition of an importan feature. If **bandskpoints** are not set
in inputs, the **BandgapWorkChain** will anyway calculate the bands following these rules:

* If a single-point calculation is requested, the kpoints path for bands is set automatically using SeeK-path.
  Please note that this choice might change the structure, as explained in the 
  `SeeK-path`_ documentation.

* If a relaxation was asked, first a siesta calculation without bands is performed to take
  care of the relaxation, then a separate single-point calculation is set up and the bands are
  calculated for a symmetry path in k-space decided by SeeK-path using the output structure of the relaxation.
  This overcomes the problem of the compatibility between bands and variable-cell relaxations.
  In fact, the final cell obtained from a relaxation, can not be known in advance, and to set
  the kpoint path without knowing the cell is generally a poor choice.
  Again note that SeeK-path might change the structure. In this second case, only the structure
  of the final single-point calculation will be changed. The changed structure is returned as
  **output_structure** port of the workchain.

If the **bandskpoints** is set by the user in inputs, no action is
taken and the behaviour follow what explained for the **SiestaBaseWorkChain**
:ref:`documentation <siesta-base-wc-inputs>`.

An additional input is present:

* **seekpath_dict** class :py:class:`Dict <aiida.orm.Dict>`, *Optional*
 
  Dictionary hosting the parametrs to pass to the ``get_explicit_kpoints_path``
  method of SeeK-path.
  The default sets ``{'reference_distance': 0.02, 'symprec': 0.0001}``,
  meaning target distance between neighboring k-points of 
  0.02 1/ang and symmetry precision parameter of 0.0001.
  Full list of the possible options and their explanation
  can be found `here`_.

Outputs
-------

* All the outputs of **SiestaBaseWorkChain** are also outputs of this 
  WorkChain, they can be explored in the relative section of the **SiestaBaseWorkChain**.

.. |br| raw:: html

    <br />
  
* **band_gap_info** :py:class:`Dict <aiida.orm.Dict>`
  
  A dictionary containing a bool (`is_insulator`) set to True if the material has a band gap,
  to False otherwise. Moreover the dictionary contains the value of the gap in `eV`.


Protocol system
---------------

The protocol system is available for this WorkChain. The ``BandgapWorkChain.inputs_generator()``
makes available all the methods explained in the :ref:`protocols documentation <how-to>`. The bands
options are still valid and they will set a `bandskpoints` input to the workchain. To avail of
the automatic generation of bands path, do not pass any ``bands_path_generator`` to ``get_filled_builder``.

.. _SeeK-path: https://seekpath.readthedocs.io/en/latest/
.. _here: https://seekpath.readthedocs.io/en/latest/module_guide/index.html#seekpath.getpaths.get_explicit_k_path
