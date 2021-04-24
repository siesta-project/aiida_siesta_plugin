"""Define warnings that can be thrown by aiida_siesta."""
from aiida.common.warnings import AiidaDeprecationWarning


class AiidaSiestaDeprecationWarning(AiidaDeprecationWarning):
    """
    Class for AiiDA siesta deprecations.
    It does *not* inherit, on purpose, from `DeprecationWarning` as
    this would be filtered out by default.
    Enabled by default, you can disable it by running in the shell::
      verdi config warnings.showdeprecations False
    """
