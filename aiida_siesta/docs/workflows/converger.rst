Converger workflow
++++++++++++++++++

Description
-----------

The **SiestaConverger** is a tool to facilitate convergence tests with Siesta.
It extends the **SiestaIterator** to accept a target quantity that is checked
after each step to evaluate whether convergence has been reached or not.
The convergence check just consists in calculating the difference in the target quantity 
between the present step and the step before and comparing it with a threshold value
passed by the user in input.
An example on the use of the **SiestaConverger** is
`/aiida_siesta/examples/workflows/example_convergence.py`.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the      
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.


.. _siesta-converger-inputs:

Inputs
------

All the **SiestaIterator** inputs are as well inputs of the **SiestaConvereger**,
they are described in the corresponding
:ref:`documentation <siesta-iterator-inputs>`.
Additional inputs are:

* **target**, class :py:class:`Str  <aiida.orm.Str>`, *Optional*

  The parameter the user wants to track in order to check if convergence has been reached.
  All the quantities returned in the **output_parameters** dictionary of the **SiestaBaseWorkChain**
  are accepted for this scope, excluding keys that don't have a `float` or `int` as a value.
  Typical values are the Kohn-Sham
  (``E_KS``), Free (``FreeE``), Band (``Ebs``), and Fermi (``E_Fermi``)
  energies, and the total spin (``stot``); however the user might also think to converge
  calculations-time related quantities.

  The `E_KS` is the default value.


.. |br| raw:: html

    <br />

* **threshold**, class :py:class:`Float <aiida.orm.Float>`, *Optional*

  The maximum difference between two consecutive steps to consider that convergence is reached.
  Default is ``Float(0.01)``.

Outputs
-------

The following outputs are returned:

* **converged** :py:class:`Bool <aiida.orm.Bool>`

  Returning `True` or `False`, whether the target has converged or not.

.. |br| raw:: html

    <br />

* **converged_target_value** :py:class:`Float <aiida.orm.Float>`

  The value of the target when the convergence has been reached. Returned only if
  the convergence is succesfull.

.. |br| raw:: html

    <br />

* **converged_parameters** :py:class:`Dict <aiida.orm.Dict>`

  The values for the parameters that was enough to achieve convergence.
  If converged is not achieved, it won't be returned.

Protocol system
---------------

The protocol system is not directly available for this WorkChain.
However inputs of the **SiestaBaseWorkChain** can be obtained in a dictionary in this way::

        inp_gen = SiestaBaseWorkChain.inputs_generator()
        inputs = inp_gen.get_inputs_dict(structure, calc_engines, protocols)

The inputs of ``get_inputs_dict`` are explained in the :ref:`protocols documentation <how-to>`.
Then the user must define at least the input **iterate_over** in order to be able to submit
the **SiestaConverger** WorkChain (if no **target** is specified, the `E_KS` is used).
