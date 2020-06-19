import numpy as np

from aiida.engine import calcfunction
from aiida.orm import Float, Str, List, Bool, Int
from aiida.plugins import DataFactory

from .iterate import SiestaIterator


@calcfunction
def generate_convergence_results(iteration_keys, variable_values, target_values, converged, converged_index):
    '''
    Generates the final output of the convergence workflows.
    '''

    convergence_results = {
        'converged': Bool(converged),
    }

    if converged:

        converged_parameters = DataFactory('dict')(dict={
            key: val for key, val in zip(iteration_keys, variable_values[converged_index.value])
        })

        convergence_results['converged_parameters'] = converged_parameters
        convergence_results['converged_target_value'] = Float(target_values[converged_index.value])

    return convergence_results


class BaseConvergencePlugin:
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
        spec.output('converged_target_value', help="The value of the target with convergence reached.")

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
                self.ctx.converged_index = below_thresh[0] + 1

            self.report(
                f'Convergence criterium: {self.inputs.threshold.value}; Last step diffs: {diffs[-len(self.ctx.last_step_processes):]}'
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

        outputs = generate_convergence_results(iteration_keys, variable_values, target_values, converged, converged_index)

        if converged:
            self.report(f'\n\nConvergence has been reached! Converged parameters: {outputs["converged_parameters"].get_dict()}\n')
        else:
            self.report('\n\nWARNING: Workchain ended without finding convergence\n ')

        self.out_many(outputs)

        super().return_results()


class SiestaConverger(BaseConvergencePlugin, SiestaIterator):
    """
    Converges a parameter in SIESTA.

    The parameter specification works exactly the same as `SiestaIterator`
    """
