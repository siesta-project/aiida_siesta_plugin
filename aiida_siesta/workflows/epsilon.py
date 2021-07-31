from aiida import orm
from aiida.engine import calcfunction
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.utils.epsilon import get_epsilon_from_eps2


@calcfunction
def get_epsilon(optical_eps2):
    """
    From an ArrayData instance, get epsilon_1(w=0)

    :param optical_eps2: the ArrayData instance with eps2 vs e

    :return epsilon: a Float instance
    """
    epsilon = get_epsilon_from_eps2(optical_eps2)

    return orm.Float(epsilon)


class EpsilonWorkChain(SiestaBaseWorkChain):
    """
    Workchain to obtain the electronic contribution to the static
    dielectric constant of a structure using Siesta.
    """

    @classmethod
    def define(cls, spec):
        """
        Add the output for the low-frequency dielectric constant and impose the input "optical" as mandatory
        """
        super().define(spec)

        spec.output('epsilon', valid_type=orm.Float, help='Low-frequency dielectric constant', required=False)
        spec.inputs["optical"].required = True

    def postprocess(self):
        """
        Calculate the low-frequency dielectric constant and return it.
        """
        super().postprocess()
        self.report('Calculating epsilon zero')
        optical_eps2 = self.outputs['optical_eps2']
        epsilon = get_epsilon(optical_eps2)
        self.out('epsilon', epsilon)
        self.report(f'EpsilonWorkChain completed. epsilon={epsilon.value}')

    @classmethod
    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
        from aiida_siesta.utils.protocols_system.input_generators import EpsilonWorkChainInputGenerator
        return EpsilonWorkChainInputGenerator(cls)
