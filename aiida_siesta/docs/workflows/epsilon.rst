Epsilon workflow
++++++++++++++++

Description
-----------

The **EpsilonWorkChain** is a simple extension of the **SietaBaseWorkChain** 
that introduces a post-processing step to obtain the low frequency dielectric
constant from the epsilon_2(omega) data.
For developers, this workchain can be taken as an example to understand how easy is to include simple
post-processes on top of the **SietaBaseWorkChain**.

Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the
`wiki <https://github.com/siesta-project/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.


Inputs
------

All the **SiestaBaseWorkChain** inputs are as well inputs of the **EpsilonWorkChain**,
therefore the system and DFT specifications (structure, parameters, etc.) are
inputted in the WorkChain using the same syntax explained in the **SiestaBaseWorkChain**
:ref:`documentation <siesta-base-wc-inputs>`.
Here we only impose a mandatory definition of the **optical** input port.

Outputs
-------

* All the outputs of **SiestaBaseWorkChain** are also outputs of this 
  WorkChain, they can be explored in the relative section of the **SiestaBaseWorkChain**.

.. |br| raw:: html

    <br />
  
* **epsilon** :py:class:`Float <aiida.orm.Float>`
  
  The low frequency (static) dielectric constant computed from the eps2(omega) data
  using Kramers-Kronig relations .


Protocol system
---------------

The protocol system is available for this WorkChain. The ``EpsilonWorkChain.inputs_generator()``
makes available all the methods explained in the :ref:`protocols documentation <how-to>`. In addition,
the **optical** input is populated, setting the optical mesh equal to the kpoints mesh of the calculation,
the "optical-broaden" to 0.5 eV and the "optical-polarization-type" to "polarized" with optical vector
of [1.0 0.0 0.0].

