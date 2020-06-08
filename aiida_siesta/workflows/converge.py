import numpy as np

from aiida.engine import WorkChain, calcfunction, while_, ToContext
from aiida.orm import Float, Str, List, KpointsData, Bool, Int
from aiida.common import AttributeDict

from .iterate import BaseIteratorWorkChain ,ParameterIterator, BasisParameterIterator,  \
                AttributeIterator, KpointComponentIterator, get_iterator_and_defaults

@calcfunction
def generate_convergence_results(variable_values, target_values, converged):
    '''
    Generates the final output of the convergence workflows.
    '''

    convergence_results = {
        'converged': Bool(converged),
    }

    if converged:
        final_value = variable_values[-2]
        dtype = {str: Str, int: Int, float: Float}[type(final_value)]

        convergence_results['final_value'] = dtype(final_value)
        convergence_results['final_target_value'] = Float(target_values[-2])

    return convergence_results

@calcfunction
def generate_new_kpoints(old_kpoints, component, val):

    mesh, offset = old_kpoints.get_kpoints_mesh()

    mesh[component.value] = val

    new_kpoints = KpointsData()
    new_kpoints.set_kpoints_mesh(mesh, offset)

    return new_kpoints

class BaseConvergenceWorkchain(BaseIteratorWorkChain):
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
        spec.input(
            "threshold",
            valid_type=(Int,Float),
            required=True,
            default=lambda: Float(0.01),
            help="The maximum difference between two consecutive steps to consider that convergence is reached"
        )

        spec.output('converged',
            help="Whether the target has converged"
        )
        spec.output('final_value',
            help="The value for the variable that was enough to achieve convergence. "
            "Otherwise, this value makes no sense."
        )
        spec.output('final_target_value',
            help="The value of the target with convergence reached."
        )
        spec.output('attempted_values')
        spec.output('target_values')

    @property
    def converged(self):
        '''
        Checks if the target property has already converged.
        '''
        target_values = self.ctx.target_values

        if len(target_values) <= 1:
            return False
        else:
            converged = abs(target_values[-2] - target_values[-1]) < self.inputs.threshold.value
            self.report(f'Convergence criterium: {self.inputs.threshold.value}; Current diff: {abs(target_values[-2] - target_values[-1])}')
            return converged

    @property
    def should_proceed(self):
        return not self.converged
    
    def process_calc(self):
        # Append the value of the target property for the last run
        results = self.ctx.calculation.outputs

        simulation_outputs = results.output_parameters.get_dict()

        self.ctx.target_values.append(simulation_outputs[self.inputs.target.value])

    def return_results(self):
        '''
        Takes care of returning the results of the workchain to the user
        '''

        converged = Bool(self.converged)
        variable_values = List(list=self.ctx.variable_values)
        target_values = List(list=self.ctx.target_values)
        
        # Return the convergence results
        self.out_many(generate_convergence_results(variable_values, target_values, converged))

        # And the 'log' of the process (this nodes have already been registered because
        # they are inputs of generate_convergence_results, which is a calcfunction)
        self.out('attempted_values', variable_values)
        self.out('target_values', target_values)

class ParameterConvergence(BaseConvergenceWorkchain, ParameterIterator):
    pass

class BasisParameterConvergence(BaseConvergenceWorkchain, BasisParameterIterator):
    pass

class AttributeConvergence(BaseConvergenceWorkchain, AttributeIterator):
    pass

class KpointComponentConvergence(BaseConvergenceWorkchain, KpointComponentIterator):
    pass

class KpointsConvergence(WorkChain):
        
    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Define the outline of the workflow, i.e. the order in which methods
        # should be executed
        spec.outline(
            cls.initialize,
            while_(cls.next_kpoint)(
                cls.find_convergence,
                cls.process_convergence_results,
            ),
            cls.return_results
        )

        # Define all the inputs that this workchain expects

        # Inputs related to the variable parameter
        spec.input(
            "order",
            valid_type=List,
            required=False,
            default=lambda: List(list=[0,1,2]),
            help="The order in which the kpoint components should be converged."
        )
        
        spec.expose_inputs(KpointComponentConvergence, exclude=('component',))
        spec.inputs._ports['pseudos'].dynamic = True
        
        # Define what are the outputs that this workchain will return
        spec.output('converged_kpoints',
            help="The kpoints with all components converged"
        )
        spec.output('final_target_value')
    
    def initialize(self):

        self.ctx.convergence_inputs = AttributeDict(self.exposed_inputs(KpointComponentConvergence))

        self.ctx.kpoints = KpointsData()
        init_value = self.ctx.convergence_inputs.init_value

        self.ctx.kpoints.set_kpoints_mesh([init_value.value]*3)

        self.ctx.left_to_converge = self.inputs.order.get_list()

    def next_kpoint(self):
        '''
        Checks if there are components left to converge.

        If there are, prepares things for next step
        '''

        if len(self.ctx.left_to_converge) == 0:
            return False

        self.ctx.component = self.ctx.left_to_converge.pop(0)

        self.report(f'Starting to converge kpoint component {self.ctx.component}')

        return True
        
    def find_convergence(self):
        '''
        Executes the corresponding convergence workflow
        '''

        # There are two inputs that are provided by this workflow
        # the rest are handled by KpointComponentConvergence
        inputs = {
            **self.ctx.convergence_inputs,
            'component': Int(self.ctx.component),
            'kpoints': self.ctx.kpoints
        }

        # Run the convergence workflow
        process = self.submit(KpointComponentConvergence, **inputs)

        # And make the workchain wait for the results
        return ToContext(last_process=process)

    def process_convergence_results(self):

        convergence_output = self.ctx.last_process.outputs

        if not convergence_output['converged']:
            raise Exception(f"Component {self.ctx.component} didn't converge")

        # If it did converge, update the kpoints
        self.report(f'Kpoint component {self.ctx.component} converged at {convergence_output["final_value"].value}')
        self.ctx.kpoints = generate_new_kpoints(self.ctx.kpoints, Int(self.ctx.component), convergence_output['final_value'])

    def return_results(self):

        self.out('converged_kpoints', self.ctx.kpoints)
        self.out('final_target_value', self.ctx.last_process.outputs['final_target_value'])

def iterator_to_converger(IteratorClass):
    
    for Converger in [ParameterConvergence, BasisParameterConvergence, AttributeConvergence]:
        if IteratorClass in Converger.__mro__:
            return Converger

def get_converger_and_defaults(iterate_over):

    Iterator, defaults = get_iterator_and_defaults(iterate_over)

    return iterator_to_converger(Iterator), defaults

def converge(over, *args, call_method='run', **kwargs):
    from aiida.engine import run, submit

    Converger, defaults = get_converger_and_defaults(over)

    execute = run if call_method is 'run' else submit

    return execute(Converger, *args, **{**defaults, **kwargs})





