SIESTA Equation Of State workflow
++++++++++++++++++++++

Description
-----------

The **EqOfStateFixedCellShape** WorkChain is a tool for the calculation of 
the equation of state of a solid.
Density Functional Theory (DFT) calculations with the SIESTA code are performed at 
7 equidistant volumes around a starting volume in order to obtain E(V) data. 
The starting volume is an optional input of the WorkChain, called `volume_per_atom`.
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
However the option relaxation `MD.ConstantVolume` (see SIESTA manual)
might be added into the parameters
`Dict <aiida.orm.Dict>` to let SIESTA relax the structure at fixed volume.
There is no point, for obvious reasons, to run this WorkChain with the 
relaxation option `MD.VariableCell`.
This WorkChain also tries to perform a Birch_Murnaghan fit
on the calculated E(V) data, following the Delta Test implementation.
If the fit fails, a warning is stored in the report of the WorkChain
(accessible through `verdi process report <PK>`), but the E(V) data for the 7 volumes 
are always returned, leading to a succesfull termination of the process.

Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1 series, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).

Inputs
------

* All the inputs of the **SiestaBaseWorkChain**, as explained 
  :ref:`here <siesta-base-wc-inputs>`.

* **volume_per_atom**,  class :py:class: `Float  <aiida.orm.Float>`, *Optional* 

  A decimal number corresponding to the volume per atom around which to 
  perform the equation of state. 

Outputs
-------

* **results_dict** :py:class:`Dict <aiida.orm.Dict>` 

  A dictionary containing a key `eos_data` containing the computed E(V) values and relative 
  units of measure.
  If the Birch-Murnaghan fit is succesfull, the key `fit_res` will report the minimum
  volume (Vo, in ang^3/atom), the minimum energy (Eo, in eV/atom), the Bulk Modulus 
  (Bo, in ev/ang^3) and its derivative respect to the presure B1.

* **equilibrium_structure** :py:class:`StructureData <aiida.orm.StructureData>`
  
  Present only if the Birch-Murnaghan fit is succesfull, it is the AiiDA structure
  at the equilibrium volume Vo.



