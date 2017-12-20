Installation
++++++++++++

It would be a good idea to create and switch to a new python virtual
environment before the installation.

Install the plugin by executing, from the top level of the plugin
directory:

::

    pip install -e .

As a pre-requisite, this will install an appropriate version of the
aiida_core python framework.

Next, run the following command

::

   reentry scan -r aiida


Development tests
+++++++++++++++++

It is possible now to run development tests located in
``aiida_siesta/tests/`` via *pytest*. The approach was originally
implemented by **Dominik Gresch** in his
`aiida\_pytest <https://github.com/greschd/aiida_pytest>`__ module.

In order to run tests, after the installation described above, run:

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
