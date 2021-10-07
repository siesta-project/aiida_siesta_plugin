import warnings
from aiida_siesta.utils.warn import AiidaSiestaDeprecationWarning
from .protocols_system.generator_absclass import InputGenerator

message = (  #pylint: disable=invalid-name
    'This module has been deprecated and will be removed in `v2.0.0`. `InputsGenerator` is now `InputGenerator` ' +
    'in `aiida_siesta.utils.protocols_sysyem.generator_absclass`.'
)

warnings.warn(message, AiidaSiestaDeprecationWarning)


class InputsGenerator(InputGenerator):  #pylint: disable=abstract-method
    """
    Mirror of new class to keep back compatibility. Since abstract class, can't be instanciated.
    """
