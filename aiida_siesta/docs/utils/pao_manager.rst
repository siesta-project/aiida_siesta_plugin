PAO manager
+++++++++++

EXPERIMENTAL FEATURE!

Description
-----------

Class to help modifications of PAO basis block. Also translates orbitals info
contained in ion files into a PAO block.
For the moment can only treat one single site at the time and can be initialized only
from an IonData instance::

        ion=load_node(pk) # pk of an IonData instance
        pao_manager = ion.get_pao_modifier()

It offers several methods to manipulate the basis specifications, for instance adding
and removing orbitals (also polarized), increase or decrease all the radii of a percentage,
manually modify single radii of orbitals.
It the returns a string that can be directly insered into the basis input of a **SiestaCalculation**
(or workchains of the package) under the `%pao-basis block` key::

        pao_manager.get_pao_block()
