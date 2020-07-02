Installation
++++++++++++

It would be a good idea to create and switch to a new python virtual
environment before the installation.

To install the plugin, simply run::

    pip install aiida-siesta

This will install the latest release of the package. Alternatively,
the source code can be downloaded from the github repository 
(https://github.com/albgar/aiida_siesta_plugin) and the insallation
can be performed by executing from the top level of the plugin
directory::

    pip install -e .

As a pre-requisite, both commands above will install an appropriate version of the
`aiida-core` python framework, if this is not already installed.
In case of a fresh install of `aiida-core`, follow the `AiiDA documentation`_
in order to configure aiida.

.. important:: In any case, do not forget to run the following commands after the installation::
                reentry scan -r aiida
                verdi daemon restart


.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/stable/
