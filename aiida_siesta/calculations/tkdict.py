import warnings
from aiida_siesta.utils.warn import AiidaSiestaDeprecationWarning
from aiida_siesta.utils.tkdict import TKDict as _TKDict
from aiida_siesta.utils.tkdict import FDFDict as _FDFDict

message = (  #pylint: disable=invalid-name
    'This module/class has been deprecated and will be removed in `v2.0.0`. ' +
    '`FDFDict` and `TKDict` should be imported from `aiida_siesta.utils.tkdict`.'
)

warnings.warn(message, AiidaSiestaDeprecationWarning)


class TKDict(_TKDict):  #pylint: disable=abstract-method
    """
    Mirror to TKDict in utiles. Do not need to put a warning since it is an abstract class
    and instanciating it will give anyway error
    """


class FDFDict(_FDFDict):
    """
    Mirror to FDFDict that rises warning when instanciated.
    """

    def __init__(self, inp_dict=None):
        warnings.warn(message, AiidaSiestaDeprecationWarning)
        super().__init__(inp_dict)
