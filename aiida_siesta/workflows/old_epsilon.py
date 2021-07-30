from aiida import orm
from aiida.engine import WorkChain, calcfunction, ToContext
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


class EpsilonWorkChain(WorkChain):
    """
    Workchain to obtain the electronic contribution to the static
    dielectric constant of a structure using Siesta.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(SiestaBaseWorkChain, exclude=('metadata', 'optical'), namespace='scf_and_relax')

        # This is required
        # The plugin will use defaults for all options
        spec.input('optical', valid_type=orm.Dict, help='Optical calculation options')

        # A prettier scheme for the user (but if we use protocols it is no better
        # than the previous one).

        # spec.input_namespace(
        #     'optical',
        #     dynamic=True,
        #     required=False,
        #     help='The inputs that will be passed to `comput_dos` method.'
        # )
        # spec.input('optical.smearing', valid_type=Float, non_db=True, required=False)
        # spec.input('optical.e_max', valid_type=Float, non_db=True, required=False)
        # spec.input('optical.e_min', valid_type=Float, non_db=True, required=False)
        # spec.input('optical.mesh', valid_type=List, non_db=True, required=False)
        # spec.input('optical.vector', valid_type=r, non_db=True, required=False)

        spec.output('optical_eps2', valid_type=orm.ArrayData, help='Array representing eps2(w)')
        spec.output('epsilon', valid_type=orm.Float, help='Low-frequency dielectric constant')

        # In case we have extra useful data
        spec.expose_outputs(SiestaBaseWorkChain, include=('output_structure', 'bands'))

        spec.outline(
            # This is a one-shot calculation, with optional relaxation, and
            # a final analysis stage that includes optical properties
            # (and maybe more analysis... see the results section)
            cls.run_siesta_wc,
            cls.run_results,
        )
        spec.exit_code(200, 'ERROR_MAIN_WC', message='The main SiestaBaseWorkChain failed')

    def run_siesta_wc(self):
        """
        Run the SiestaBaseWorkChain, which might include relaxation
        """

        inputs = self.exposed_inputs(SiestaBaseWorkChain, namespace='scf_and_relax')

        # Add optical block, which was not passed from the workchain
        inputs['optical'] = self.inputs.optical

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report(f'Launched SiestaBaseWorkChain<{running.pk}> (scf [+relax] + optical).')

        return ToContext(workchain_base=running)

    def run_results(self):

        if not self.ctx.workchain_base.is_finished_ok:
            return self.exit_codes.ERROR_FINAL_WC

        outps = self.ctx.workchain_base.outputs

        optical_eps2 = outps['optical_eps2']

        epsilon = get_epsilon(optical_eps2)

        self.out('optical_eps2', optical_eps2)
        self.out('epsilon', epsilon)

        # Get other results
        # (this might be the seed of a design for an "omnibus" analysis workchain)

        if 'output_structure' in outps:
            self.out('output_structure', outps['output_structure'])
        if 'bands' in outps:
            self.out('bands', outps['bands'])

        self.report(f'EpsilonWorkChain completed. epsilon={epsilon.value}')


#    @classmethod
#    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
#        from aiida_siesta.utils.inputs_generators import BaseWorkChainInputsGenerator
#        return BaseWorkChainInputsGenerator(cls)
