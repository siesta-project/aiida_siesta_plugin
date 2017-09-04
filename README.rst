AiiDA Siesta plugins and workflows
==================================

This repository contains the files that implement the AiiDA Siesta
plugins and a set of related AiiDA workflows.

Documentation can be found in:

http://aiida-siesta-plugin.readthedocs.io

Installation and usage guidelines
---------------------------------

Installation
~~~~~~~~~~~~

The distribution package of the plugin is available on PyPI.
To install it directly from the index run:

::

       pip install --process-dependency-links aiida-siesta

However, at this development stage ``aiida\_siesta`` uses a custom version of ``aiida\_core``, to work around a few issues relevant to the Siesta plugin. Thus, using package from PyPI is not recommended at the moment, since the installation relies on deprecated ``--process-dependency-links`` flag.

For a more suitable development setup (including custom ``aiida\_core`` installation) see the section below.


Local environment setup
^^^^^^^^^^^^^^^^^^^^^^^

Use your own installation of the database, and your own computers and
codes.

Assuming that you are
at the top-level of the Siesta AiiDA plugin hierarchy, do:

::

       git clone https://github.com/vdikan/aiida_core
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

Docker-based usage and automated testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Docker framework is being updated to the new plugin architecture.
