Installation
++++++++++++

It would be a good idea to create and switch to a new python virtual
environment before the installation.

The latest release of the package can be obtained simply with::

    pip install aiida-siesta

In this case, make sure to refer to the appropriate documentation part ("stable", not "latest").

Because the package is under development, in order to enjoy the most recent features
one can clone the github repository
(https://github.com/albgar/aiida_siesta_plugin) and install
from the top level of the plugin directory with::

    pip install -e .

As a pre-requisite, both commands above will install an appropriate version of the
``aiida-core`` python framework, if this is not already installed.
In case of a fresh install of ``aiida-core``, follow the `AiiDA documentation`_
in order to configure aiida.

.. important:: In any case, do not forget to run the following commands after the 
   installation::
                
        reentry scan -r aiida
        verdi daemon restart


For developers
--------------

This plugin is open-source and contributions are welcomed. Before starting the development, the following steps
are suggested:

* After cloning from github, install with ``pip install .[dev]``. This will download all the tools for testing.
* Install `pre-commit`_ hooks. This will "force" to follow some python standards we require. In fact, the hooks will impede 
  to commit unless the required standards are met.
* Make sure to run all the tests (simply ``pytest test/`` from the main folder of the package) to make sure the contribution is not 
  breaking any part of the code. Ideally, write tests for the new part implemented.

.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/stable/
.. _pre-commit: https://pre-commit.com/#install
