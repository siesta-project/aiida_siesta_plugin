Bandgap workflow
++++++++++++++++

Description
-----------

The **BandgapWorkChain** is an extension of the **SietaBaseWorkChain** 
that introduces a simple post-process with the scope to return the metallic or
insulating nature of the material and, possibly, the band gap.
The purpose of this WorkChain is mostly educational, showing how easy is
to introduce pre-processes or post-processes in the WorkChain logic.
The class **BandgapWorkChain** is, in fact, a subclass of the **SietaBaseWorkChain**
that just overrides few methods and introduces the
additional output **band_gap_info**.

To calculate the gap, this workchain makes use of a tool distributed in aiida-core,
the method ``find_bandgap`` hosted in ``aiida.orm.nodes.data.array.bands``.

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
The only difference is that the **bandskpoints** are now a mandatory input and the WorkChain
will rise an error if they are not present.

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
makes available all the methods explained in the :ref:`protocols documentation <how-to>`, the
only difference is that here is mandatory to pass ``bands_path_generator`` to ``get_filled_builder`` and
not optional like for the **SietaBaseWorkChain** inputs generator.
