.. AiiDA Siesta Plugin documentation master file, created by
   sphinx-quickstart on Mon Jun 26 11:33:57 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the AiiDA-Siesta documentation!
++++++++++++++++++++++++++++++++++++++++++++

..

The aiida-siesta python package interfaces the SIESTA DFT code
(https://siesta-project.org/siesta/) with the AiiDA framework
(http://www.aiida.net).  The package contains: plugins for SIESTA
itself and for other utility programs, new data structures, and basic
workflows. It is distributed under the MIT license and available from
(https://github.com/siesta-project/aiida_siesta_plugin).
If you use this package, please cite J. Chem. Phys. **152**, 204108 (2020)
(https://doi.org/10.1063/5.0005077).


Acknowledgments:
----------------

The Siesta input plugin was originally developed by Victor M. Garcia-Suarez.

Alberto Garcia further improved the Siesta input plugin and wrote the parser for Siesta and the STM plugin.

Emanuele Bosoni contributed the band-structure support for the Siesta plugin.

Vladimir Dikan and Alberto Garcia developed the workflows and
refined the architecture of the package.

Vladimir Dikan and Emanuele Bosoni ported the plugin and the base workflow to AiiDA 1.0.
Alberto Garcia futher refined the 1.0-compatible functionality.

From November 2019 to May 2022, Emanuele Bosoni was in charge of the code's development and maintenance,
under the supervision of Alberto Garcia.

Pol Febrer contributed the SiestaIterator and SiestaConverger workflows, including the underline
abstract classes system.

This work is supported by the `MaX
European Centre of Excellence <http://www.max-centre.eu/>`_ funded by
the Horizon 2020 INFRAEDI-2018-1 program, Grant No. 824143, and by the
`INTERSECT <https://intersect-project.eu/>`_  (Interoperable material-to-device simulation box for
disruptive electronics) project, funded by Horizon 2020 under grant
agreement No 814487, as well as by the Spanish MINECO (projects
FIS2012-37549-C05-05 and FIS2015-64886-C5-4-P)

We thank the AiiDA team, who are also supported by the `MARVEL National Centre for Competency in Research <http://nccr-marvel.ch>`_
funded by the `Swiss National Science Foundation <http://www.snf.ch/en>`_

.. figure:: miscellaneous/logos/MaX.png
    :alt: MaX
    :target: http://www.max-centre.eu/
.. figure:: miscellaneous/logos/INTERSECT.png
    :alt: INTERSECT
    :target: http://intersect-project.eu/
.. figure:: miscellaneous/logos/MINECO-AEI.png
    :alt: MINECO-AEI
    :target: http://www.mineco.gob.es/
.. figure:: miscellaneous/logos/MARVEL.png
    :alt: MARVEL
    :target: http://nccr-marvel.ch

Contents:
---------

Installation
============

.. toctree::
   :maxdepth: 4

   installation

Calculation plugins
===================

.. toctree::
   :maxdepth: 4

   plugins/index


Utilities
=========

.. toctree::
   :maxdepth: 4

   utils/index


Workflows
=========

.. toctree::
   :maxdepth: 4

   workflows/index
..

Tutorials
=========

.. toctree::
   :maxdepth: 4

   tutorials/index
