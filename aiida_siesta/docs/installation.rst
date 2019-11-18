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

.. important:: Next, do not forget to run the following commands

::

   reentry scan -r aiida
   verdi daemon restart

