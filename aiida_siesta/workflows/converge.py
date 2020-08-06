import numpy as np

from aiida.engine import calcfunction
from aiida.orm import Float, Str, List, Bool, Int
from aiida.plugins import DataFactory

from .iterate import SiestaIterator, InputIterator, BaseIteratorWorkChain


@calcfunction
def generate_convergence_results(iteration_keys, variable_values, target_values, converged, converged_index):
    '''
    Generates the final output of the convergence workflows.
    '''

    convergence_results = {
        'converged': Bool(converged),
    }

    if converged:
        #pylint: disable=unnecessary-comprehension
        converged_parameters = DataFactory('dict')(
            dict={key: val for key, val in zip(iteration_keys, variable_values[converged_index.value])}
        )

        convergence_results['converged_parameters'] = converged_parameters
        convergence_results['converged_target_value'] = Float(target_values[converged_index.value])

    return convergence_results


class BaseConvergencePlugin(InputIterator):
    '''
    Plugin to add to an Iterator workchain to convert it to a convergence workflow.

    To use it, just build a class that inherits from this plugin and the iterator
    that you want to use.

    Examples
    -----------
    class MyConverger(BaseConvergencePlugin, MyIterator):
        pass

    This class will now iterate as MyIterator while checking for convergence.

    This class just checks between the difference in a target output between two
    consecutive steps and compares it to a threshold. To implement a different
    convergence algorithm, just overwrite the `converged` property.
    '''

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # The outline is defined by Iterator workflow, and most of the inputs too.
        # Threshold is an input that is specific to a convergence workflow.
        # Inputs related to the parameter to converge
        spec.input(
            "target",
            valid_type=Str,
            required=False,
            default=lambda: Str('E_KS'),
            help="The parameter that you want to track."
        )
        spec.input(
            "threshold",
            valid_type=(Int, Float),
            required=True,
            default=lambda: Float(0.01),
            help="The maximum difference between two consecutive steps to consider that convergence is reached"
        )

        spec.output('converged', help="Whether the target has converged")
        spec.output(
            'converged_parameters',
            required=False,
            help="The values for the parameters that was enough to achieve convergence. "
            "If converged is not achieved, it won't be returned",
        )
        spec.output('converged_target_value', required=False, help="The value of the target with convergence reached.")

    def initialize(self):
        """
        Initialize the list with values of the target output.
        """
        super().initialize()

        self.ctx.target_values = []

    @property
    def converged(self):
        '''
        Checks if the target property has already converged.
        '''
        target_values = self.ctx.target_values

        if len(target_values) <= 1:
            converged = False
        else:
            diffs = abs(np.diff(target_values))
            below_thresh = np.where(diffs < self.inputs.threshold.value)[0]

            converged = len(below_thresh) > 0
            if converged:
                self.ctx.converged_index = below_thresh[0]

            self.report(
                'Convergence criterium: '
                '{0}; Last step diffs: {1}'.format(
                    self.inputs.threshold.value, diffs[-len(self.ctx.last_step_processes):]
                )
            )

        return converged

    def should_proceed(self):
        return not self.converged

    def analyze_process(self, process_node):
        """
        Takes the output from the process and stores the value of the target.
        """
        # Append the value of the target property for the last run
        results = process_node.outputs

        simulation_outputs = results.output_parameters.get_dict()

        self.ctx.target_values.append(simulation_outputs[self.inputs.target.value])

        super().analyze_process(process_node)

    def return_results(self):
        '''
        Takes care of returning the results of the convergence to the user
        '''

        converged = Bool(self.converged)
        converged_index = Int(getattr(self.ctx, 'converged_index', -1))
        iteration_keys = List(list=list(self.ctx.iteration_keys))
        variable_values = List(list=self.ctx.variable_values)
        target_values = List(list=self.ctx.target_values)

        outputs = generate_convergence_results(
            iteration_keys, variable_values, target_values, converged, converged_index
        )

        if converged:
            self.report(
                '\n\nConvergence has been reached! Converged parameters:'
                f'{outputs["converged_parameters"].get_dict()}\n'
            )
        else:
            self.report('\n\nWARNING: Workchain ended without finding convergence\n ')

        self.out_many(outputs)

        super().return_results()


class SiestaConverger(BaseConvergencePlugin, SiestaIterator):
    """
    Converges a parameter in SIESTA.

    The parameter specification works exactly the same as `SiestaIterator`

    Examples
    ------------

    See `examples/workflows/example_convergence.py`.
    """


class SequentialConverger(InputIterator):
    '''
    Launches several convergence workchains sequentially.

    At each step, it incorporates the already converged parameters from the previous one.

    THIS CLASS CAN NOT BE USED DIRECTLY, you need to subclass it and specify a cls._process_class,
    which should be the converger that you want to use.
    '''

    _process_class = None
    _expose_inputs_kwargs = {'exclude': ("iterate_over",), "namespace": 'converger_inputs'}
    _reuse_inputs = True

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # In this case, iterate_over is going to be a list where each item is a dict as accepted
        # BaseIterator's regular iterate_over
        spec.inputs.ports['iterate_over'].valid_type = List

        # In principle we should not allow to provide a batch size and we should force a value of 1
        # But I don't know how to do it in a clean way.
        spec.inputs.ports['batch_size'].default = lambda: Int(1)

        spec.output(
            'converged_parameters',
            required=False,
            help="The values for the parameters that was enough to achieve convergence. "
            "If converged is not achieved, it won't be returned",
        )
        spec.output('converged_target_value', help="The value of the target with convergence reached.")

    @classmethod
    def _iterate_input_serializer(cls, iterate_over):
        """
        Parses the "iterate_over" key of the workchain.

        Parameters
        -----------
        iterate_over: list of dict or aiida Dict

            For each dictionary, each key is the name of a parameter we want to iterate
            over (str) and each value is a list with all the values to iterate over for
            that parameter.

        Returns
        -----------
        aiida Dict
            the parsed input.
        """

        if isinstance(iterate_over, list):
            parsed_list = []
            for step in iterate_over:
                parsed_list.append(BaseIteratorWorkChain._iterate_input_serializer(step))
            return cls._values_list_serializer(parsed_list)

        return iterate_over

    def initialize(self):
        super().initialize()

        self.ctx.already_converged = {}
        self.ctx._input_keys = {}
        self.ctx._parsing_funcs = {}

    def _get_iterator(self):
        """
        Builds the iterator that will generate values.
        """

        iterate_over = self.inputs.iterate_over.get_list()

        self.ctx.iteration_keys = ("iterate_over",)

        return zip(iterate_over)

    def analyze_process(self, process_node):

        if process_node.outputs.converged:

            # Get the parameters that have resulted in convergence
            new_converged = process_node.outputs.converged_parameters.get_dict()

            # Store them in what is going to be the final output
            self.ctx.already_converged.update(new_converged)

            # And now we are going to modify the inputs for the next step according
            # to what has converged, so that the convergence of the next parameter already incorporates
            # what we have found.
            for key, val in new_converged.items():

                # Get the last inputs
                inputs = self.ctx.last_inputs

                # Modify them accordingly. Note that since we are using cls._reuse_inputs = True,
                # what we are modifying here will be used in the next iteration
                self._params_lookup = self._process_class._params_lookup
                self._process_class._parse_key(self, key)
                val = self._process_class._parse_val(self, val, key, inputs)

                key = self._process_class._attr_from_key(self, key)

                self._process_class.add_inputs(self, key, val, inputs)

                self._params_lookup = self.__class__._params_lookup

        self.ctx.last_target_value = process_node.outputs.converged_target_value

    def return_results(self):

        self.out('converged_parameters', DataFactory('dict')(dict=self.ctx.already_converged).store())
        self.out('converged_target_value', self.ctx.last_target_value)


class SiestaSequentialConverger(SequentialConverger):
    _process_class = SiestaConverger
