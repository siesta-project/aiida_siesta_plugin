AiiDA Siesta plugins and workflows
==================================

A plugin with `Siesta DFT code <https://departments.icmab.es/leem/siesta/>`_
interface to the `AiiDA system <http://www.aiida.net/>`_.

Documentation can be found in:

http://aiida-siesta-plugin.readthedocs.io

Installation and usage guidelines
---------------------------------

Installation
~~~~~~~~~~~~

``aiida-siesta`` is compatible with latest release of AiiDA.

To install it directly from PyPI run:

::

       pip install aiida-siesta


Development guide
~~~~~~~~~~~~~~~~~

For development and testing purposes, create and activate python ``virtualenv`` environment instance, then follow the steps below:


Local environment setup
^^^^^^^^^^^^^^^^^^^^^^^

::

       git clone https://github.com/aiidateam/aiida_core
       cd aiida_core
       git checkout develop

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
^^^^^^^^^^^^^^^^^

It is possible now to run development tests located in
``aiida_siesta/tests/`` via *pytest*.
The approach we use for unit- and integrity testing is implemented
by **Dominik Gresch** in his `aiida\_pytest <https://github.com/greschd/aiida_pytest>`__ module.

In order to run tests, after you install core and plugin as described above, run:

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

Docker-based usage and automated testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Docker framework is being updated to the new plugin architecture.
