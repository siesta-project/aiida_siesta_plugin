Equation Of State workflow
++++++++++++++++++++++++++

Description
-----------

The **EqOfStateFixedCellShape** WorkChain is a tool for the calculation of 
the equation of state of a solid.
Density Functional Theory (DFT) calculations with the SIESTA code are performed at 
7 equidistant volumes around a starting volume in order to obtain the energy (E)
versus volume (V) data.
The starting volume is an optional input of the WorkChain, called **volume_per_atom**.
If the latter is not specified, the input structure volume is use as starting volume.
The WorchChain ensure robustness in the convergence of each SIESTA calculation thanks to 
the fact that each DFT run is submitted through the **SiestaBaseWorkChain**,
that automatically manages some common failures (lack of
electronic-structure or geometry relaxation convergence, termination due to
walltime restrictions, etc).
All the **SiestaBaseWorkChain** inputs are as well inputs of the **EqOfStateFixedCellShape**,
therefore the system and DFT specifications (structure, parameters, etc.) are
inputted in the WorkChain using the same syntax explained in the **SiestaBaseWorkChain**
:ref:`documentation <siesta-base-wc-inputs>`.
As the name of the class suggest, the **EqOfStateFixedCellShape** is designed to
obtain the E(V) curve under the restriction of fixed cell shape.
This means that no algorithm for stress minimization is implemented in the WorkChain.
However the option relaxation ``MD.ConstantVolume`` (see SIESTA manual)
might be added into the parameters
dictionary to let SIESTA to relax the structure at fixed volume.
There is no point, for obvious reasons, to run this WorkChain with the 
relaxation option ``MD.VariableCell``.
This WorkChain also tries to perform a Birch_Murnaghan fit
on the calculated E(V) data, following the `DeltaProject`_ implementation.
If the fit fails, a warning is stored in the report of the WorkChain
(accessible through ``verdi process report <PK>``), but the E(V) data for the 7 volumes 
are always returned, leading to a succesfull termination of the process.

Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the      
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.



Inputs
------

* All the inputs of the **SiestaBaseWorkChain**, as explained 
  :ref:`here <siesta-base-wc-inputs>`.

.. |br| raw:: html

    <br />

* **volume_per_atom**, class :py:class:`Float <aiida.orm.Float>`, *Optional* 

  A decimal number corresponding to the volume per atom around which to 
  perform the equation of state. 

.. |br| raw:: html

    <br />

* **batch_size**, class :py:class:`Int <aiida.orm.Int>`, *Optional*

  Number of volumes to run at the same time. By default, it is set to one,
  therefore one volume at the time is submitted


Outputs
-------

* **results_dict** :py:class:`Dict <aiida.orm.Dict>` 

  A dictionary containing a key `eos_data` that collects the computed E(V) values and relative 
  units of measure.
  If the Birch-Murnaghan fit is succesfull, also the key `fit_res` will be present in this disctionary.
  It reports the following values extracted from the fit: the equilibrium
  volume (Vo, in `Angstom^3/atom`), the minimum energy (Eo, in `eV/atom`), the Bulk Modulus 
  (Bo, in `ev/Angstrom^3`) and its derivative respect to the presure B1.  
 
.. |br| raw:: html

    <br />
  
* **equilibrium_structure** :py:class:`StructureData <aiida.orm.StructureData>`
  
  Present only if the Birch-Murnaghan fit is succesfull, it is the AiiDA structure
  at the equilibrium volume Vo.



Protocol system
---------------

The protocol system is available for this WorkChain. The ``EqOfStateFixedCellShape.inputs_generator()``
makes available all the methods explained in the :ref:`protocols documentation <how-to>`, the
only difference is that the relaxation type "variable-cell" is not available.


.. _DeltaProject: https://github.com/molmod/DeltaCodesDFT/blob/master/eosfit.py
