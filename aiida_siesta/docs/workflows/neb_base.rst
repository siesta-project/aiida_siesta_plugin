NEB Base workflow
++++++++++++++++

Description
-----------
The **SiestaBaseNEBWorkChain** is the core building block for the creation of workflows that
enable the search of the Minimum Energy Pathway (MEP) connecting two local minima of the potential
energy surface through the Nudge Elastic Band (NEB) method.
In particular, this workchain performs NEB MEP optimizations starting from a guessed path
and exploiting the LUA functionality of SIESTA.
This workchain is very useful for the investigation of reaction paths and energy barriers. For instance,
it can be used to study the energetic barrier for interstitial diffusion of an impurity in
an host structure. For some concrete examples, look at the `aiida-siesta-barrier` project.

An example on the use of the **SiestaBaseNEBWorkChain** is in
`/aiida_siesta/examples/workflows/example_neb_ghost.py`.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, 4.1-b3 of the 4.1 series and the MaX-1.0 release, which
can be found in the development platform
(https://gitlab.com/siesta-project/siesta).
For more up to date info on compatibility, please check the
`wiki <https://github.com/siesta-project/aiida_siesta_plugin/wiki/Supported-siesta-versions>`_.

.. _siesta-base-neb-wc-inputs:

Inputs
------

Many of the **SiestaBaseWorkChain** inputs are as well inputs of the **SiestaBaseNEBWorkChain**,
therefore the system and DFT specifications (structure, parameters, etc.) can be defined as
input in the WorkChain using the same syntax explained in the **SiestaBaseWorkChain**
:ref:`documentation <siesta-base-wc-inputs>`.
The only exceptions are the **structure** and the **lua** namespace
that are not explicit inputs for this workchain.
In fact, more than one single structure is required by NEB method and they are passed
through the dedicated input **starting_path**. The **lua** inputs are mostly defined
internally except the lua script that is now named **neb_script**. A more detailed
description of the two new inputs follows:

* **starting_path**, class :py:class:`TrajectoryData <aiida.orm.TrajectoryData>`, *Mandatory*

  A set of structures collected in a
  :py:class:`TrajectoryData <aiida.orm.TrajectoryData>` object. Each structure correspond to
  an image for the NEB method. The object must have the kinds of the structure as
  attributes.

 .. |br| raw:: html

    <br />

* **neb_script**, class :py:class:`SingleFileData <aiida.orm.SingleFileData>`, *Mandatory*

  A lua script that controls the NEB calculation. An example can be seen in
  `/aiida_siesta/examples/fixtures/lua_scripts/neb.lua`.

.. note::  The use of LUA scripts also requires the user to pass to aiida the environmental
   variable that indicates where the flos library is. More info :ref:`here <submission-siesta-calc>`.

Outputs
------

* **neb_output_package**, class :py:class:`TrajectoryData <aiida.orm.TrajectoryData>`, *Mandatory*

  A :py:class:`TrajectoryData <aiida.orm.TrajectoryData>` object with the final structures after
  the NEB optimization and the energy of each one of them. Moreover the reaction barrier and
 other useful info are reported as attributes of the node.



Protocol system
---------------
No protocol system is in place for this workchain.
