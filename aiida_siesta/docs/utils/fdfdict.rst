FDF dictionary
++++++++++++++

Description
-----------

The FDFDict class represents data from a .fdf-file (the standard input of the siesta
code). It behaves like a normal python
dictionary, but with translation rules that follow the standards of the Flexible Data Format (FDF).
The FDF format was developed inside the siesta package in order to facilitate the
creation of the input file of siesta. Among other features, it substitute strings in favour of
default values.
In particular it drops dashes/dots/colons and imposes lowercase.
The FDFDict class accepts in input a python dictionary and applies the same
rules to the "keys" of the dictionary.
An example::

         from aiida_siesta.calculations.tkdict import FDFDict
         inp_dict = {"ThisKey": 3,"a-no-ther": 4,"t.h.i.r.d" : 5}
         f = FDFDict(inp_dict)
         print(f.keys())

returns ``dict_keys(['thiskey', 'another', 'third'])``.

When two keys in the same dictionary will become the same string after translation, the last
definition will remain::
        
         from aiida_siesta.calculations.tkdict import FDFDict
         inp_dict = {"w":3,"e":4,"w--":5}
         f = FDFDict(inp_dict)
         print(f.get_dict())

returns ``{'w': 5, 'e': 4}``.

The method ``get_dict`` returns the translated dictionary, but the class keeps record also of
the last unstraslated key for each key.
This can be seen just printing ``f``. The method
``get_untranslated_dict`` returns the dictionary with the last unstranslated keys as keys.
Therefore in our example, the ``get_untranslated_dict`` returns ``{'w--': 5, 'e': 4}``.

Getter and setter are implemented to get and set the value automatically for each equivalent
key. ``f["w"]``, ``f["w---"]`` will return the same value. The call ``f["w---"] = 3`` will reset
the value of key ``"w"``, also changing the "last untranslated key" to ``"w---"``.

Many more methods are available in the FDFDict class. They can be explored from the source code
(``aiida_siesta.calculations.tkdict``).
It is a useful tool for the development of new CalcJobs and WorkChains.


