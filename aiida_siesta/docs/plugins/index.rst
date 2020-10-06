Calculations
^^^^^^^^^^^^

This section contains the documentation for the calculations plugins 
distributed in ``aiida-siesta``.
They are the fundamental blocks that enable to run through AiiDA some executable of the Siesta package,
meaning the siesta code itself and some post-processing tools.
For each calculation, we explain the inputs selection, the submission command and the returned outputs.
From the AiiDA prospective, we describe here the functionalities of both the CalcJob class and 
the associated parser. 

.. toctree::
   :maxdepth: 2

   siesta
   stm
