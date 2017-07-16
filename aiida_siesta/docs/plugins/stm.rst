STM  plugin
++++++++++++++++++++++

Description
-----------

A plugin for Util/plstm


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1
series, which can be found in the development platform (http://launchpad.net/siesta/).

Inputs
------

* **parameters**, class :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`

A dictionary with a few parameters to specify the mode of calculation
and the height or isovalue at which to process the LDOS::

    {
      "mode": "constant-height",
      "z": "5.8"     # In Ang
    }

* **parent_folder**, class
  :py:class:`RemoteData <aiida.orm.data.RemoteData>`
      
The parent folder of a previous Siesta calculation in which the LDOS
file was generated.

Outputs
-------


* **stm_array** :py:class: `ArrayData <aiida.orm.data.array.ArrayData>` 

A 2D array holding the section or topography information.

* **output_parameters** :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>` 
  (accessed by ``calculation.res``)

At this point only parser information is returned.


Errors
------

Errors during the parsing stage are reported in the log of the calculation (accessible 
with the ``verdi calculation logshow`` command). 


