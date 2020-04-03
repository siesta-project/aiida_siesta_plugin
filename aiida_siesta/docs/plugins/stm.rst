STM  plugin
++++++++++++++++++++++

Description
-----------

A plugin for Util/plstm of the Siesta distribution, a tool to simulate STM images.


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1
series, which can be found in the development platform (https://gitlab.com/siesta-project/siesta)

Inputs
------

Some examples are referenced in the following list. They are located in the folder aiida_siesta/examples/plugins/stm.

* **code**, class :py:class:`Code <aiida.orm.Code>`, *Mandatory*

  A code object linked to a plstm executable.
  If you setup the code ``plstm1`` on machine ``kelvin`` following the `aiida guidelines`_,
  then the code is selected in this way::

        codename = 'plstm1@kelvin'
        from aiida.orm import Code
        code = Code.get_from_string(codename)


* **mode**, class :py:class:`Str <aiida.orm.Str>`, *Mandatory*

  Allowed values are `constant-height` or `constant-current`, corresponding to the two
  operation modes of the STM that are supported by the `plstm` code.
  Examples for both modes are presented in the example folder.

.. |br| raw:: html

    <br />

* **value**, class :py:class:`Float <aiida.orm.Float>`, *Mandatory*

  The value of height or current at which the user wants to simulate the 
  STM. The height must be expressed in Ang, the current in e/bohr**3.

.. |br| raw:: html

    <br />


* **ldos_folder**, class :py:class:`RemoteData <aiida.orm.RemoteData>`, *Mandatory*
      
  The parent folder of a previous Siesta calculation in which the .LDOS
  file was generated. To have more information on how to produce the .LDOS file,
  one can refer to the example aiida_siesta/examples/plugins/siesta/example_ldos.
  Please note that the `ldos_folder` must be on the same machine on which the STM analysis
  is performed. In other words, the input **code** must be installed on the same machine 
  where the `ldos_folder` resides. This is a limitation of AiiDA that can not copy
  between different computers, but it is also required by `plstm` itself, as the .LDOS
  file is produced in an unformatted way.

.. |br| raw:: html

    <br />

* **spin_option**, class :py:class:`Str <aiida.orm.Str>`, *Optional*

  Input port that allows the selection of the spin options offered by `plstm`. It follows the same
  syntax of the code. The value "q" selects a total charge analysis. The value "s" selects the 
  total spin magnitude analyisis (only available if the parent Siesta calculation is spin polarized).
  Finally, the values "x", "y" or "z" indicate a separate analysis of one the three spin components
  (only available if the parent Siesta calculation is performed with non-collinear options).
  If the port is not specified the default "q" option is activated.

.. |br| raw:: html

    <br />

* **settings**, class :py:class:`Str <aiida.orm.Str>`, *Optional*

  A port `settings` is available to activate some advanced features. For instance the modification
  of the command line instructions and the addition of files to retreave. For more info,
  the corresponding section of the Standard Siesta Plugin can be seen :ref:`here <siesta-advanced-features>`.


Submitting the calculation
--------------------------

The submission of any CalcJob of AiiDA always follows the same schema. Therefore,
to understand how to submit a STM calculation, it is sufficient to follow the explanation
of the corresponding section of the Standard Siesta Plugin.
The only change is to import the correct plugin::

        from aiida_siesta.calculations.stm import STMCalculation
        builder = STMCalculation.get_builder()
and, of course, to define the correct inputs allowed by `STMCalculation` (previous 
section).


Outputs
-------

* **stm_array** :py:class:`ArrayData <aiida.orm.ArrayData>` 

  A collection of three 2D arrays (`grid_X`, `grid_Y`, `STM`) holding the section or
  topography information. They follow the `meshgrid` convention in
  Numpy. A contour plot can be generated with the `get_stm_image.py`
  script in the repository of examples.

.. |br| raw:: html

    <br />

* **output_parameters** :py:class:`Dict <aiida.orm.Dict>` 

  At this point, it constains only the parser information and the name of the 
  retrieved file where the STM info were stored.


Errors
------

Errors during the parsing stage are reported in the log of the calculation (accessible 
with the ``verdi process report`` command). 


