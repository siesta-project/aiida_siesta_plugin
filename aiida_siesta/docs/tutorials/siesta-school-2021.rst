.. _Siesta school 2021 Homepage:

2021, Virtual event
===================


+-----------------+----------------------------------------------------------------------------+
| Related resources                                                                            |
+=================+============================================================================+
| Virtual Machine | `Siesta Mobile 0.2.0`_                                                     |
+-----------------+----------------------------------------------------------------------------+
| python packages | `aiida-core 1.6.1`_, `aiida-siesta 1.2.0`_,                                |
+-----------------+----------------------------------------------------------------------------+
| codes           | `Siesta Max-1.2.0-1`_                                                      |
+-----------------+----------------------------------------------------------------------------+

.. _Siesta Mobile 0.2.0: https://drive.google.com/drive/u/2/folders/14V50YRuJfW1jxdWkQzZPnTx0TIa10ftX
.. _aiida-core 1.6.1: https://pypi.org/project/aiida-core/1.6.1
.. _aiida-siesta 1.2.0: https://pypi.org/project/aiida-siesta/1.2.0
.. _Siesta Max-1.2.0-1: https://gitlab.com/siesta-project/siesta/-/wikis/Guide-to-Siesta-versions

These are the notes of the tutorial delivered to the CECAM school "First-principles simulations of materials with SIESTA"
running virtually from 28th of June to 2nd of July 2021.
The tutorial was carried on using the Siesta Mobile Virtual Machine, however reference to the relevant aiida documentation
is reported at the beginning. If you are running on Siesta Mobile, jump :ref:`here <start>` 
Tutors: Emanuele Bosoni

Installation and setup (ONLY IF NOT IN SIESTA MOBILE)
-----------------------------------------------------

Installation is through ``pip`` after moving to a new virual environment (we use ``virtualenvwrapper``, but any alternative is valid, only
make sure to select a python version 3.6 or above). We call the virtual environment ``tutorial``.

.. code-block:: python

   mkvirtualenv tutorial
   workon tutorial
   pip install aiida==1.6.1
   pip install aiida-siesta==1.2.0

Follow the instructions in the AiiDA
`documentation <https://aiida.readthedocs.io/projects/aiida-core/en/v1.6.1/intro/get_started.html>`_.
to set up aiida.

Computer and code setup (ONLY IF NOT IN SIESTA MOBILE)
------------------------------------------------------

Follow the `instructions <https://aiida.readthedocs.io/projects/aiida-core/en/v1.6.1/howto/run_codes.html>`_
to set up your computer and a siesta executable.

.. _start:

Creating a pseudo family
------------------------

In the Siesta Mobile, activate the virtual environment ``workon siesta_school``.
Otherwise activate the environment you created before.

Check the status of aiida typing ``verdi status``. Check the codes installed with
``verdi code list``.

Before starting to play with `aiida-siesta`, it can be useful to learn how to
set up a pseudopotential family containing a collection of pseudos from `PseudoDojo <http://www.pseudo-dojo.org/>`_.
Just do::

     aiida-pseudo install pseudo-dojo -v 0.4 -x PBE -r SR -p standard -f psml

This will install version 0.4 PBE scalar relativistic and standard accuracy in psml form under the name
"PseudoDojo/0.4/PBE/SR/standard/psml".

If you have a FOLDER containing psf pseudopotentials, you can create a family with::

     aiida-pseudo install family /PATH/TO/FOLDER/ FAM_NAME -P pseudo.psf


Submit a single siesta calculation
----------------------------------

Open the file :download:`example_bands.py <data_siesta_school_2021/example_bands.py>` and explore the setting up of the various inputs.
Run the script with (if not in Siesta Mobile use the approppriate code name you set up)::

        runaiida example_bands.py --dont-send siesta-school--Max-1.3.0-1@localhost

The option ``--dont-send`` has been added in order to activate the "dry_run" option that every aiida process has.
This option allows to create all the inputs of the calculation, but do not submit it.
You can explore in the folder ``submit_test`` how AiiDA prepared all the inputs of a siesta calculation for you.

Now run::

        runaiida example_bands.py --send siesta-school--Max-1.3.0-1@localhost

AiiDA took charge of your script, created the inputs and submitted the calculation. Look at the state of the
process with the command ``verdi process show <pk>`` as suggested in the shell. The ``<pk>`` number
uniquely identify your calculation and it will be used later on.

In few seconds the calculation is finished. You will relized that when ``verdi process show <pk>`` 
shows "Finished" status and reports the oututs. We explore the outputs. This can be done from command line,
for instance::

        verdi data array show <PK_forces_and_stress>

        
however it is worth exploring the shell provided by AiiDA::

        verdi shell

Inside the shell::

        l=load_node(<PK_calculation>)

and explore all the methods making use of tab completion. For instance::

        l.outputs.bands.export(path="Si_bands", fileformat="gnuplot", y_max_lim=10)

The command above creates a file that can be plot with gnuplot in order to visualize the bands. Open a new shell and
type::

        gnuplot --persist Si_bands

Take the chance to explore in the ``verdi shell`` some methods and attributes of data types associated to 
the inputs and outputs of a ``SiestaCalculation``. Use tab complition of ``l.inputs``, ``l.outputs``, ``l.attributes``, ..

The submission script can be modified very easily in order to run a ``SiestaBaseWorkChain`` instead of a
``SiestaCalculation``. Look at the commented part of the ``example_bands.py`` script in order to understand the differences.
The "dry_run" option is not available for the ``SiestaBaseWorkChain``.
A ``SiestaBaseWorkChain`` automatically takes care
of fixing some `common errors of a siesta calculation
<https://aiida-siesta-plugin.readthedocs.io/en/latest/workflows/base.html#error-handling>`_,
therefore it adds robustness in running siesta calculations.

Protocols
---------
Go back to the ``verdi shell`` and look at the following::

        from aiida_siesta.workflows.base import SiestaBaseWorkChain
        inp_gen=SiestaBaseWorkChain.inputs_generator()

You just imported the inputs generator for the ``SiestaBaseWorkChain``. We can explore its functionality::

        inp_gen.get_protocol_names()
        inp_gen.get_spins()

And many more... Use tab complition to explore them. These methods allows you
to understand which options you can pass to ``get_filled_builder``, as will be explained in a second.

The main feature of the input generator is the possibility to obtain a ``builder`` (a tool that helps you build 
the inputs for the specific process) that is ready to be submitted::

        struct = l.inputs.structure
        calc_engines = {
            'siesta': {
                'code': "siesta-school--Max-1.3.0-1@localhost",
                'options': {'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1},"max_wallclock_seconds": 3600}
                }
            }
        builder = inp_gen.get_filled_builder(struct,calc_engines,"standard_psml")

The ``calc_engines`` is a dictionary with fixed keys, whose aim is to pass the computational resourses
for the calculation.

Explore the ``builder``::

        builder.parameters.attributes
        builder.basis.attributes
        ...

We can add spin polarization to the calculation with::

        builder = inp_gen.get_filled_builder(struct,calc_engines,"standard_psml",spin="polarized")

Try again ``builder.parameters.attributes``, what are the differences compared to before?

We can run the builder straight away::

        from aiida.engine import run
        run(builder)

Wait the calculation to finish and explore the results if you wish.

We are now going to create our own protocol.
Close the shell and look at the file :download:`my_protocols_registry.yaml <data_siesta_school_2021/my_protocols_registry.yaml>`. 
It contains a new set of inputs and the psml pseudos. 
This file can be modified at will and its content will become a new protocol. Simply look at the folder where
you are ``pwd`` and attach the file to the correct environment variable, like that:: 

        export AIIDA_SIESTA_PROTOCOLS="path_discovered_with_pwd/my_protocols_registry.yaml"

taking care of passing the correct absolute path where you have ``my_protocols_registry.yaml``.

Now open the shell and::

        from aiida_siesta.workflows.base import SiestaBaseWorkChain
        inp_gen=SiestaBaseWorkChain.inputs_generator()
        inp_gen.get_protocol_names()

The new protocol is on the list and we can use it to run a calculation::

        l=load_node(<PK_calculation>)
        struct = l.inputs.structure
        calc_engines = {
            'siesta': {
                'code': "siesta-school--Max-1.3.0-1@localhost",
                'options': {'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1},"max_wallclock_seconds": 3600}
                }
            }
        builder = inp_gen.get_filled_builder(struct,calc_engines,"my_protocol")

        from aiida.engine import run
        run(builder)

The command ``run`` send the calcualation in the shell in interactive mode (does not submit to the builder
as ``submit`` would do).
Our set up will occupy the shell for a minute or so and at the end it will return the outputs of the calculation.


Run a convergence workflow
--------------------------

It's quite easy to run a convergence workflow using `aiida-siesta`.

In a verdi you can do::

        from aiida_siesta.workflows.converge import SiestaSequentialConverger
        from aiida.engine import run

        calc_node=load_node(<PK_calculation>)

        run(SiestaSequentialConverger,

                iterate_over=[
                        {"kpoints_0": [4,5,6,7,8,9,10,11,12,13,14,15]},
                        {"kpoints_1": [4,5,6,7,8,9,10,11,12,13,14,15]}
                ],

                converger_inputs={
                        'code':load_code('siesta-school--Max-1.3.0-1@localhost'),
                        'pseudo_family': Str('siesta-school--Max-1.3.0-1@localhost')
                        'structure': calc_node.inputs.structure,
                        'parameters': Dict(),
                        'options': {'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1},"max_wallclock_seconds": 3600}
                        'batch_size': Int(3)
                }

        )

to converge your structure's kpoints (first and second components), running three simulations at a time.
The ``<PK_calculation>`` is the integer we got at the beginning of the tutorial.

Want to know more??
--------------------------

In general about Aiida (create your workflos and so on)? `Aiida tutorials <https://aiida-tutorials.readthedocs.io/en/latest/>`_

On aiida siesta? `docs <https://docs.siesta-project.org/projects/aiida-siesta/en/latest/index.html>`_

Ask me: ebosoni@icmab.es
