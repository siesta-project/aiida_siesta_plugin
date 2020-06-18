from aiida.engine import WorkChain, calcfunction, while_, ToContext
from aiida.orm import Float, Str, List, KpointsData, Bool, Int
from aiida.plugins import DataFactory
from aiida.common import AttributeDict

from .iterate import SiestaIterator


@calcfunction
def generate_convergence_results(iteration_keys, variable_values, target_values, converged):
    '''
    Generates the final output of the convergence workflows.
    '''

    convergence_results = {
        'converged': Bool(converged),
    }

    if converged:

        final_value = DataFactory('dict')(dict={
            key: val for key, val in zip(iteration_keys, variable_values[-2])
        })

        convergence_results['final_value'] = final_value
        convergence_results['final_target_value'] = Float(target_values[-2])

    return convergence_results


class BaseConvergencePlugin:
    '''
    General workflow that finds the converged value for a parameter.

    There is no specificity in what should be converged or which parameter
    should be varied.
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
            'final_value',
            help="The value for the variable that was enough to achieve convergence. "
            "Otherwise, this value makes no sense."
        )
        spec.output('final_target_value', help="The value of the target with convergence reached.")

    def initialize(self):
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
            converged = abs(target_values[-2] - target_values[-1]) < self.inputs.threshold.value
            self.report(
                f'Convergence criterium: {self.inputs.threshold.value}; Current diff: {abs(target_values[-2] - target_values[-1])}'
            )

        return converged

    @property
    def should_proceed(self):
        return not self.converged

    def analyze_process(self, process_node):
        # Append the value of the target property for the last run
        results = process_node.outputs

        simulation_outputs = results.output_parameters.get_dict()

        self.ctx.target_values.append(simulation_outputs[self.inputs.target.value])

    def return_results(self):
        '''
        Takes care of returning the results of the Plugin to the user
        '''

        converged = Bool(self.converged)
        iteration_keys = List(list=list(self.ctx.iteration_keys))
        variable_values = List(list=self.ctx.variable_values)
        target_values = List(list=self.ctx.target_values)

        # Return the convergence results
        outputs = generate_convergence_results(iteration_keys, variable_values, target_values, converged) 
        self.out_many(outputs)


class SiestaConverger(BaseConvergencePlugin, SiestaIterator):
    pass
