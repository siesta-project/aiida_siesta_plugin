Basis optimization
++++++++++++++++++++

Description
-----------

The AiiDA-siesta package offers three workchains to help users in selecting the optimal
basis set for a given system:

1) The **SimplexBasisOptimization** that finds the minimum basis enthalpy varying a set of input variables (typically 
cutoff radii of orbitals) using the Nelder–Mead (simplex / amoeba) method.
2) The **TwoStepsBasisOpt** that performs a two level optimization, running simplex iterations followed by periodic 
restarts with new simplex hyper-tetrahedra of progressively smaller sizes. This replicates more closely the
optimization util of the siesta distribution.
3) The **BasisOptimizationWorkChain** that performs a full optimization testing first basis cardinality and then applying
the **SimplexBasisOptimization** to optimize all the radius of the orbitals.

#An example on the use of the **SiestaConverger** is
#`/aiida_siesta/examples/workflows/example_seq_converger.py`


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the      
`wiki <https://github.com/siesta-project/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.



SimplexBasisOptimization
------------------------

The Nelder–Mead method (commonly known as simplex or amoeba method) is 
a numerical method to find the minimum of a function with N variables.
A simplex is a special polytope of N+1 vertices in N dimensions 
(for instance a triangle in 2D, a tetrahedron in 3D and so forth). The
Nelder–Mead methods in N dimensions starts from a set of N+1 test points arranged as a simplex. 
The value of the function is calculated at each test point and these values are then used
in order to find a new test point and 
to replace one of the old test points with the new one. Repeating the procedure, the technique progresses
until all the N+1 test points produce a value within a threshold. When this happens the minimum has been reached.
In the context of the basis optimization the function is the basis enthalpy and the variable
are the parameters defining the basis (typically cutoff radii of the basis orbitals).

The simplex optimization is performed taking advantage of the ``aiida-optimize`` package (link), in particular the
Nelder–Mead method implemented in that package.

Inputs
******

Inputs are organized in two namespaces and are described in the following:

* **siesta_base**, input namespace, *Mandatory*
  
  Accepts all the inputs of a **SiestaBaseWorkChain** (listed `here <siesta-base-wc-inputs>`) with the only mandatory modification
  to include in the "basis" input some variables to optimize. The variable must be defined using a dollar and
  a string. An example::

       basis = Dict(dict={
        '%block pao-basis': "\nSi   2\n n=3   0   2\n 4.99376      $sz2 \n n=3   1   2 P 1\n 6.2538      $pz2 \n%endblock pao-basis"
        }) 

  An upper and lower value must be set for each variable and optionally one or more starting points (see next point in this list). 
  Please note that variables are typically defined for orbitals radii in the ``pao-basis`` block,
  but one can also create variables for "higher level" keywords like the ``energy-shift`` or ``split-norm``.

.. |br| raw:: html

    <br />


* **simplex**, input namespace, *Mandatory*
  
  Here all the inputs for the simplex method can be defined. They are listed in the next lines.

.. |br| raw:: html

    <br />

* **simplex.variables_dict** class :py:class:`Dict <aiida.orm.Dict>`, *Mandatory*

  A dictionary containing all the info about the variables that are modified in order to find the minimum
  basis enthalpy. An example related to the basis block above::

        variables_dict = Dict(dict={
            "sz2":[2.0,4.8,3.0],
            "pz2":[2.0,6.0,3.0]
            })

  The variables names must be defined here as keys of the dictionary and must correspond to the
  strings defined in the ``basis`` input, but removing the dollar symbol.
  The list associated to each string defines in this order: 1) The lower limit for the variable,
  2) The upper limit, 3) the starting value to construct the simplex hyper-tetrahedron.
  The up and down limit of the variables are used in such way: if the algorithm attempts
  the calculation of the function out of range, a huge value for the function is returned.
  The starting value is going to be the point from which the simplex hyper-tetrahedron is constructed.
  In particular, the first test point is directly formed by the specified starting points (in the example above is [3.0,3.0]).
  The other N test points are obtained substituing one component with ``num + range *  simplex_inps.initial_step_fraction``,
  where ``num`` is the defined starting point, ``range`` is the upper - lower limit and ``simplex_inps.initial_step_fraction`` 
  is a number between 0 and 1 defined in the next point of this list.
  Supposing ``simplex_inps.initial_step_fraction = 0.2``, in out example, the other two test points are [3.0,3.8] and 
  [3.56,3.0].

  When 3) is not defined, it is chosen randomly between the boundaries, but it is always suggested
  to set it since it will be used to construct the 
  Alternatively to 3), N+1 values can be entered and this would correspond to define explicitly all the components of
  all the simplex initial points. 

.. |br| raw:: html

    <br />


* **simplex.initial_step_fraction** class :py:class:`Float <aiida.orm.Float>`, *Optional*

  A fractional increment to be used in the construction of the starting simplex hyper-tetrahedron.
  See point above for more details. Default at ``Float(0.4)``. It is ignored if all the components
  af all the test points are set in the point above.

.. |br| raw:: html

    <br />


* **simplex.max_iters** class :py:class:`Int <aiida.orm.Int>`, *Optional*

  The maximum iterations for the Nelder–Mead algorithm. Please note that an iteration step usually involves more then one new
  test point. So the points tested at the end will be way more than the ``max_iters``.
  Once the ``simplex.max_iters`` is reached, the workchain stops returning the best simplex so far, even if the
  threshold convergence has not been reached.
  Default is ``Int(40)``.

.. |br| raw:: html

    <br />


* **simplex.tolerance_function**  class :py:class:`Float <aiida.orm.Float>`, *Optional*

  The tolerance accepted to define the optimization converged. If the values of the functions for all 
  points in the simplex are all within the ``simplex.tolerance_function``, the optimization is considered concluded.
  The default is ``Float(0.01)``.
  Please note that the choice of this parameter must be related to the variance of the output function.
  For the moment only the basis enthalpy is supported as function, therefore the default is reasonable, but future
  extensions might require to think to have a fractional tolerance.


Outputs
*******

The following outputs are returned:

* **last_simplex**  class :py:class:`List <aiida.orm.List>`

  The output containing the values of the last simplex. Always returned, even if the optimization does not reached the
  required tolerance. It is a list of lists. The first element of the list is always the best choice of the parameters
  obtained by the optimization so far.

.. |br| raw:: html

    <br />

* **optimal_process_input** class :py:class:`List <aiida.orm.List>`

  This output contains the optimal set of parameters obtained after optimization. This corresponds to the first entry of
  the list return by the **last_simplex**, however it is returned only if the optimization succeed.

.. |br| raw:: html

    <br />


* **optimal_process_output** class :py:class:`Float <aiida.orm.Float>`

  The value of the function for the optimal set of parameters obtained with the optimization.
  Returned only if the optimization succeed.

.. |br| raw:: html

    <br />


* **optimal_process_uuid** class :py:class:`List <aiida.orm.List>`

  The uuid of the **SiestaBaseWorkChain** that has the **optimal_process_input** as variables and that
  returned the **optimal_process_output**. Returned only if the optimization succeed.

It is important to note that the optimization is entirely an AiiDA process, therefore the provenance of all calculation called is preserved.
We can have a look at the attempted variables values and the obtained basis entalpy in this simple way. In the verdi shell::

        node=load_node(<PK>)  #PK of your SimplexBasisOptimization
        for wc in node.called[0].called: 
             print(wc.inputs.the_values.get_list(),wc.outputs.ene.value) 

And many more info can be extracted from the inputs and outputs of each run ``wc``. These ``wc`` are **SiestaBaseWorkChain**
wrapped into a thin layer that attach to each calculation the information needed by the optimizer.


TwoStepsBasisOpt
----------------

This workchain uses the **SimplexBasisOptimization**, but it adds a step in the optimization,
which consists in restarting the simplex with a subsequently smaller **simplex.initial_step_fraction**.
This is implemented in the original simplex optimization code that can be found
in the Util of the SIESTA package. There the fractional step is called "lambda" and we will follow the same
notation here.

Inputs
*******

All the inputs of **SimplexBasisOptimization** are inputs of this workchain except the **simplex.initial_step_fraction**.
This include the way to specify the optimization variables in the ``siesta_base.basis`` input.
This workchain adds a further called **macrostep**. This allows:

* **macrostep.initial_lambda** class :py:class:`Float <aiida.orm.Float>`

  The value of lambda to be used as **simplex.initial_step_fraction** in the first iteration.
  Default ``Float(0.4)``,


.. |br| raw:: html

    <br />

* **macrostep.lambda_scaling_factor** class :py:class:`Float <aiida.orm.Float>`
  
  The rate at which lambda decreases between from a macrostep to the other.
  Default ``Float(0.5)``


.. |br| raw:: html

    <br />

* **macrostep.minimum_lambda** class :py:class:`Float <aiida.orm.Float>`
  
  When this value for lambda is reached, the macrostep iteration stops. Default ``Float(0.01)``.



Outputs
*******

Same outputs of **SimplexBasisOptimization**.


BasisOptimizationWorkChain
--------------------------

This workchain manages entirely the optimization of the basis sets for a SIESTA calculation.
It first run calculations with different basis sizes (using the "PAO-BasisSize" option of SIESTA)
and gets the size that gives minimum basis enthalpy.

NOTE: This does not include yet the possibility to test different basis sizes for different species.

It then automatically sets up a **SimplexBasisOptimization** where the radii of all orbitals are considered variables.
The first zetas are allowed between 3 and 12 Bohrs, the subsequent zetas are defined as ratios of the first zetas
and the ratio varies between 0.2 to 0.9.

Inputs
*******

All the inputs of **SimplexBasisOptimization** are inputs of this workchain except the **simplex.variables_dict**
and the **siesta_base.basis**.
Few more inputs are allowed:

* **basis_sizes** class :py:class:`List <aiida.orm.List>` *Optional*

  The list of basis sizes to try out. Default ``List(list=["DZ", "DZP", "TZ"])``.

.. |br| raw:: html

    <br />
    
* **non_perturbative_pol** class :py:class:`Bool <aiida.orm.Bool>` *Optional*

  In order to activate the "PAO-non-perturbative-polarization-orbitals" logical keyword of SIESTA, that,
  if enabeled, promotes any polarization shells to the status of explicit shells, thus using the
  standard generation options. This is applied already in the part where the basis sizes are tested.
  Default ``Bool(False)``
  This is a fairly important feature since the "normal" (perturbative) polarization schema of SIESTA
  does not allow control on the radii of the polarized orbitals. In setting **non_perturbative_pol** to True,
  all the orbital are listed as standard orbitals in the PAO block and therefore their radii can be optimized
  like any other orbital.

Outputs
*******

Only one output is produced:

* **optimal_basis_block** class :py:class:`Dict <aiida.orm.Dict>`

  Returning the optimal pao block, meaning the one that gives minimum basis enthalpy.



Protocol system
---------------

The protocol system is not directly available for this WorkChain.
However inputs of the **SiestaBaseWorkChain** can be obtained in a dictionary in this way::

        inp_gen = SiestaBaseWorkChain.inputs_generator()
        inputs = inp_gen.get_inputs_dict(structure, calc_engines, protocols)

The inputs of ``get_inputs_dict`` are explained in the :ref:`protocols documentation <how-to>`.
Then the user can place these ``inputs`` in the **siesta_base** namespace.
