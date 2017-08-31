Installation
++++++++++++

You need to work on a specially patched copy of aiida\_core that works
around a few issues relevant to the Siesta plugin. Assuming that you
are at the top-level of the Siesta AiiDA plugin hierarchy, do:

::

       git clone https://github.com/vdikan/aiida_core
       cd aiida_core
       git checkout v0.9_siesta_patch

Install aiida with the '-e' option (and optionally install the
dependencies for extra utilities such as pymatgen, cif, etc, and the
generation of docs)

::

       cd aiida_core
       pip install -e .
       pip install -e .[docs,atomic_tools]

The first time you will need to configure AiiDA as explained in the
manual.

Install the plugin by executing, from the top level of the plugin
directory:

::

    pip install -e .

Development tests
+++++++++++++++++

It is possible now to run development tests located in
``aiida_siesta/tests/`` via *pytest*. The approach was originally
implemented by **Dominik Gresch** in his
`aiida\_pytest <https://github.com/greschd/aiida_pytest>`__ module.

In order to run tests, after you install core and plugin as described
above, run:

::

    pip install -r test_reqirements.txt

which will install necessary versions of dependencies.

Then ``cd aiida_siesta/tests`` and invoke:

::

    ./all.sh

to run all tests, or:

::

    ./run.sh <test_module_filenames>

to run specific tests.
