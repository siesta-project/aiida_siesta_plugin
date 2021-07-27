import warnings
from aiida_siesta.utils.warn import AiidaSiestaDeprecationWarning
import aiida_siesta.utils.protocols_system.input_generators as inp_gen

message = (  #pylint: disable=invalid-name
    'This module/class has been deprecated and will be removed in `v2.0.0`. ' +
    'the input generators are now in `aiida_siesta.utils.protocols_system.input_generators`.'
)

warnings.warn(message, AiidaSiestaDeprecationWarning)


class SiestaCalculationInputsGenerator(inp_gen.SiestaCalculationInputGenerator):
    """
    Mirror class that rises warning.
    """

    def __init__(self, workchain_class):
        warnings.warn(
            message + " `SiestaCalculationInputsGenerator` changed its name to `SiestaCalculationInputGenerator`",
            AiidaSiestaDeprecationWarning
        )
        super().__init__(workchain_class)


class BaseWorkChainInputsGenerator(inp_gen.BaseWorkChainInputGenerator):
    """
    Mirror class that rises warning.
    """

    def __init__(self, workchain_class):
        warnings.warn(
            message + " `BaseWorkChainInputsGenerator` changed its name to `BaseWorkChainInputGenerator`",
            AiidaSiestaDeprecationWarning
        )
        super().__init__(workchain_class)


class EosWorkChainInputsGenerator(inp_gen.EosWorkChainInputGenerator):
    """
    Mirror class that rises warning.
    """

    def __init__(self, workchain_class):
        warnings.warn(
            message + " `EosWorkChainInputsGenerator` changed its name to `EosWorkChainInputGenerator`",
            AiidaSiestaDeprecationWarning
        )
        super().__init__(workchain_class)


class StmWorkChainInputsGenerator(inp_gen.StmWorkChainInputGenerator):
    """
    Mirror class that rises warning.
    """

    def __init__(self, workchain_class):
        warnings.warn(
            message + " `StmWorkChainInputsGenerator` changed its name to `StmWorkChainInputGenerator`",
            AiidaSiestaDeprecationWarning
        )
        super().__init__(workchain_class)
