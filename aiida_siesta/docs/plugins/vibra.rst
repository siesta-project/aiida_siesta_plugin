STM  plugin
++++++++++++++++++++++

Description
-----------

A plugin for Util/Vibra


Supported Siesta versions
-------------------------

At least 4.0.1 of the 4.0 series, and 4.1-b3 of the 4.1
series, which can be found in the development platform (http://launchpad.net/siesta/).


Inputs
------

* **parameters**, class :py:class:`ParameterData <aiida.orm.data.parameter.ParameterData>`

A dictionary with the supercell parameters, the atomic displacement and
a flag to calculate or not the eigenvectors::

    {
      "supercell_1": 1, 
      "supercell_2": 1, 
      "supercell_3": 1,
      "atomicdispl": "0.021 Ang",
      "eigenvectors": True,
    }

* **bandskpoints**, class :py:class:`KpointsData <aiida.orm.data.parameter.KpointsData>`

Two lists of variables, one with the k-points (an integer which specifies
the number of k-points followed by three reals that specify their position),
and another withe the k-points labels::

        kpp = [(1,1.,1.,1.),
               (15,0.,0.5,0.5),
               (25,0.,0.,0.),
               (20,0.5,0.5,0.5),
               (20,0.,0.5,0.5),
               (15,0.25,0.5,0.75),
               (20,0.5,0.5,0.5)]
        lpp = [[0,'\Gamma'],
               [1,'X'],
               [2,'\Gamma'],
               [3,'L'],
               [4,'X'],
               [5,'W'],
               [6,'L']]

* **parent_folder**, class
  :py:class:`RemoteData <aiida.orm.data.RemoteData>`
      
The parent folder of a previous Siesta calculation in which the FC
file was generated, **OR**

* **singlefile**, class
  :py:class:`SingefileData <aiida.orm.data.SinglefileData>`
      
The local folder which contains the FC file from a previous Siesta calculation.


Outputs
-------

Two files, one with the phonon frequencies (aiida.bands) and the other
with the phonon eigenvectors (aiida.vectors).


Errors
------

Errors during the parsing stage are reported in the log of the calculation
(accessible with the ``verdi calculation logshow`` command). 


