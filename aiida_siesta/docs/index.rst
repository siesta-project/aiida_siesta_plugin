.. AiiDA Siesta Plugin documentation master file, created by
   sphinx-quickstart on Mon Jun 26 11:33:57 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the AiiDA-Siesta documentation!
++++++++++++++++++++++++++++++++++++++++++++

..

The aiida-siesta python package interfaces the SIESTA DFT code
(http://www.icmab.es/siesta) with the AiiDA framework
(http://www.aiida.net).  The package contains: plugins for SIESTA
itself and for other utility programs, new data structures, and basic
workflows.  It is distributed under the MIT license and available from
(https://github.com/albgar/aiida_siesta_plugin).

Acknowledgments:
----------------

The Siesta input plugin was originally developed by Victor
M. Garcia-Suarez.

Alberto Garcia further improved the Siesta input plugin and wrote the
parser for Siesta and the STM plugin.

Emanuele Bosoni contributed the band-structure support for the Siesta
plugin.

Vladimir Dikan and Alberto Garcia developed the workflows and
refined the architecture of the package.

We acknowledge partial support from the Spanish Research Agency (projects
FIS2012-37549-C05-05, FIS2015-64886-C5-4-P and PGC2018-096955-B-C44) and  by the `MaX 
European Centre of Excellence <http://www.max-centre.eu/>`_ funded by the Horizon 2020 EINFRA-5 program,
Grant No. 676598.

We thank the AiiDA team, who are also supported by the [MARVEL National Centre for Competency in Research](<http://nccr-marvel.ch>)
funded by the `Swiss National Science Foundation <http://www.snf.ch/en>`_

.. figure:: miscellaneous/logos/MINECO-AEI.png
    :alt: MINECO-AEI
.. figure:: miscellaneous/logos/MaX.png
    :alt: MaX
.. figure:: miscellaneous/logos/MARVEL.png
    :alt: MARVEL

Contents:
---------

Installation
============

.. toctree::
   :maxdepth: 2

   installation

SIESTA plugins
=================

.. toctree::
   :maxdepth: 2

   plugins/siesta
   plugins/stm

SIESTA Workflows
======================

.. toctree::
   :maxdepth: 4

   workflows/base
   workflows/eos
   workflows/bands
   workflows/stm
..
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

