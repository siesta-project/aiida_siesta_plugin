STM  plugin
++++++++++++++++++++++

Description
-----------

A plugin for Util/plstm


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1
series, which can be found in the development platform (https://gitlab.com/siesta-project/siesta)

Inputs
------

* **parameters**, class :py:class:`Dict <aiida.orm.Dict>`

A dictionary with a few parameters to specify the mode of calculation
and the height or isovalue at which to process the LDOS::

    {
      "z": "5.8"     # In Ang
    }

(The `mode of calculation` is hard-wired to `constant-height` for now)

* **ldos_folder**, class
  :py:class:`RemoteData <aiida.orm.RemoteData>`
      
The parent folder of a previous Siesta calculation in which the LDOS
file was generated.

Outputs
-------


* **stm_array** :py:class:`ArrayData <aiida.orm.ArrayData>` 

A collection of three 2D arrays (`X`, `Y`, `Z`) holding the section or
topography information. They follow the `meshgrid` convention in
Numpy. A contour plot can be generated with the `get_stm_image.py`
script in the repository of examples.

* **output_parameters** :py:class:`Dict <aiida.orm.Dict>` 

At this point only parser information is returned.


Errors
------

Errors during the parsing stage are reported in the log of the calculation (accessible 
with the ``verdi process report`` command). 


