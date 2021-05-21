import warnings
from aiida_siesta.utils.warn import AiidaSiestaDeprecationWarning
from .protocols_system.protocols import ProtocolManager as _ProtocolManager

message = (  #pylint: disable=invalid-name
    'This module/class has been deprecated and will be removed in `v2.0.0`. ' +
    '`ProtocolManager` should be imported from `aiida_siesta.utils.protocols_system.protocols`.'
)

warnings.warn(message, AiidaSiestaDeprecationWarning)


class ProtocolManager(_ProtocolManager):
    """
    Mirror class to print warning
    """

    def __init__(self):
        warnings.warn(message, AiidaSiestaDeprecationWarning)
        super().__init__()
