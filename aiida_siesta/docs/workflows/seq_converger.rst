Sequential Converger workflow
+++++++++++++++++++++++++++++

Description
-----------

The **SiestaSequentialConverger** is an iterator that sequentially runs **SiestaConvergers**.
Once the convergence over a parameter is reached, the converged value is used for the
following convergence test (on a new parameter).
An example on the use of the **SiestaConverger** is
`/aiida_siesta/examples/workflows/example_seq_converger.py`


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the      
`wiki <https://github.com/albgar/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.



Inputs
------

Two are the required inputs:

* **converger_inputs**, class :py:class:`dict`, *Mandatory*

  A dictionary containing all the inputs required by the **SiestaConverger**, except the 
  **iterate_over** port. The explanations of the converger inputs can be examined
  `here <siesta-converger-inputs>`. Please note that the normal inputs of a **SiestaBaseWorkChain**
  process (structure, parameters, basis, code, ...) must be included as well in this dictionary.

  The same default values as **SiestaConverger** apply if some ports are not specified here.


.. |br| raw:: html

    <br />

* **iterate_over**, class :py:class:`list`, *Mandatory*

  There is a specific port for the quantities to iterate over and now the accepted value for
  this port is a `list`, not a dictionary like it was for the **SiestaConverger** or **SiestaIterator**.
  In fact, now the user should indicate a list of parameters that he/she wants to converge
  sequentially.
  A practical example::

        iterate_over=[
           {
            'kpoints_0': [4,10,12,14,16,18,20],
            'kpoints_1': [4,10,12,14,16,18,20],
            'kpoints_2': [4,10,12,14,16,18,20],
           },
           {
            'meshcutoff': ["500 Ry", "600 Ry", "700 Ry", "800 Ry", "900 Ry"],
           },
           {
            'pao-energyshift': ["0.02 Ry", "0.015 Ry", "0.01 Ry", "0.005 Ry", "0.001 Ry"]
           }
        ]


  With this specification, we signal that we want to converge first the kpoints by increasing all components
  at the same time (assuming "zip" is selected as 'iterate_mode' in the **converger_inputs** dictionary),
  then the 'meshcutoff' and finally the 'energy shift'. The converged kpoints will be used for the convergence
  of 'meshcutoff', the converged kpoints and 'meshcutoff' will be used for the convergence process of 'energy shift'.

  Note that one can converge the same parameters again if wanted,
  for instance set up different rounds for kpoints convergence.

.. warning:: If one of the parameters does not converge, no action is taken and
   the following convergence step is performed using the inputs specified in **converger_inputs**,
   not using the last attempted value in the previous convergence.
   For instance, in the example above, if the ``meshcutoff`` does not converged at ``900 Ry``,
   the ``pao-energyshift`` convergence will be done using the inputs parameters
   specified in the ``parameters`` of ``converger_inputs``, not including
   ``meshcutoff = "900 Ry"``.



Outputs
-------

The following outputs are returned:


* **converged_target_value** :py:class:`Dict <aiida.orm.Dict>`

  The value of the target when the convergence has been reached. Returned only if
  at least one of the sequential convergences has been completed succesfull.

.. |br| raw:: html

    <br />

* **converged_parameters** :py:class:`Dict <aiida.orm.Dict>`

  The values for the parameters that was enough to achieve convergence.
  If converged is not achieved, it will be an empty dictionary.

.. |br| raw:: html

    <br />

* **unconverged_parameters** :py:class:`List <aiida.orm.List>`

  If one or more parameters fail to converge, we list them
  in this output.


Protocol system
---------------

The protocol system is not directly available for this WorkChain.
However inputs of the **SiestaBaseWorkChain** can be obtained in a dictionary in this way::

        inp_gen = SiestaBaseWorkChain.inputs_generator()
        inputs = inp_gen.get_inputs_dict(structure, calc_engines, protocols)

The inputs of ``get_inputs_dict`` are explained in the :ref:`protocols documentation <how-to>`.
Then the user can place these ``inputs`` in the **converger_inputs** dictionary (together with the other
**SiestaConverger** inputs specifications). The input **iterate_over** is also required
in order to be able to submit the **SiestaSequentialConverger** WorkChain and it must be set manually.
