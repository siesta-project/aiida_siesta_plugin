Iterator workflow
+++++++++++++++++

Description
-----------

The **SiestaIterator** is a tool to facilitate the submission of several Siesta Calculations
in an automatic way. It allows the iteration over Siesta parameters
and, more in general, over inputs of a **SiestaBaseWorkChain**.
An example on the use of the **SiestaConverger** is
`/aiida_siesta/examples/workflows/example_iterate.py`.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the      
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.


.. _siesta-iterator-inputs:

Inputs
------

All the **SiestaBaseWorkChain** inputs are as well inputs of the **SiestaIterator**,
therefore the system and DFT specifications (structure, parameters, etc.) are
inputted in the WorkChain using the same syntax explained in the **SiestaBaseWorkChain**
:ref:`documentation <siesta-base-wc-inputs>`.
The additional inputs are:

* **iterate_over**, class :py:class:`Dict  <aiida.orm.Dict>`, *Mandatory*

  A dictionary where each key is the name of a parameter we want to iterate
  over (:py:class:`str <str>`) and each value is a :py:class:`list <list>` with all the values to iterate over for
  the corresponding key.  
  Accepted keys are:

  * Name of the input ports of the **SiestaBaseWorkChain**. Meaning all the names listed
    :ref:`here <siesta-base-wc-inputs>`.
    In this case, the corresponding values list must contains the list of :py:class:`Data <aiida.orm.Data>` 
    nodes (stored or unstored) accepted by the key. Examples are::

        code1 = load_code("SiestaHere@localhost")
        code2 = load_code("SiestaThere@remotemachine")
        iterate_over = {"code" : [code1,code2]}

        struct1 = StructureData(ase=ase_struct_1)
        struct2 = StructureData(ase=ase_struct_2)
        iterate_over = {"structure" : [struct1,struct2]}

  * Name of accepted Siesta input keywords (for instance ``mesh-cutoff``, ``pao-energy-shift``, etc ...).
    In this case, the corresponding values list must contains the list of values directly, meaning
    :py:class:`str <str>`, :py:class:`float <float>`, :py:class:`int <int>` or :py:class:`bool <bool>` 
    python types. Examples are::

        iterate_over = {"spin" : ["polarized", "spin-orbit"]}

    .. warning:: In order to guarantee full flexibility, no check on the Siesta parameters is performed. If you pass as key something not recognized by Siesta, the SiestaIterator will include it in the `parameters` input and run the calculation with no warning issued. Because Siesta will not understand the keyword, it will ignore it, resulting in a series of identical calculations.
    
  The `iterate_over` is a dictionary because it is possible to iterate over several keywords at
  the same time. Something of this kind::

        struct1 = StructureData(ase=ase_struct_1)
        struct2 = StructureData(ase=ase_struct_2)
        iterate_over = {"structure" : [struct1,struct2], "spin" : ["polarized", "spin-orbit"]}

  is perfectly acceptable and the way the algorithm handle with these multiple iterations is decided
  by the **SiestaIterator** input explained next in this list.

.. |br| raw:: html

    <br />

* **iterate_mode**, class :py:class:`Str <aiida.orm.Str>`, *Optional*

  Indicates the way the parameters should be iterated. Currently allowed values are
  'zip' (zips all the parameters together, this imposes that all keys should
  have the same number of values in the list!) and 'product' (performs a cartesian product of the 
  parameters, meaning that all possible combinations of parameters and values are explored).

  The option 'zip' is the default one.

.. |br| raw:: html

    <br />

* **batch_size**, class :py:class:`Int <aiida.orm.Int>`, *Optional*

  The maximum number of simulations that should run at the same time.
  You can set this to a very large number if you want that all simulations run in
  one single batch. As default, only one single calculation at the time is submitted.


Outputs
-------

This WorkChain does not generate any output! It is, in fact, a tool to help the
submission of multiple calculations and keep them all connected and easy accessible
through the main workchain node, but it does not have any precise scope.
AiiDA provides a powerful `querying system`_ to explore all the results of the submitted calculations
and a tool to `organize the data`_.


Protocol system
---------------

The protocol system is not directly available for this WorkChain.
However inputs of the **SiestaBaseWorkChain** can be obtained in a dictionary in this way::

        inp_gen = SiestaBaseWorkChain.inputs_generator()
        inputs = inp_gen.get_inputs_dict(structure, calc_engines, protocols)

The inputs of ``get_inputs_dict`` are explained in the :ref:`protocols documentation <how-to>`.
Then the user must define at least the input **iterate_over** in order to be able to submit
the **SiestaIterator** WorkChain.

.. _querying system: https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/data.html#finding-and-querying-for-data
.. _organize the data: https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/data.html#organizing-data
