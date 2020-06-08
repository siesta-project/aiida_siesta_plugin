import itertools
from functools import partial
import numpy as np

from aiida.plugins import DataFactory
from aiida.engine import WorkChain, calcfunction, while_, ToContext
from aiida.orm import Float, Str, List, KpointsData, Bool, Int
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.common import AttributeDict

from ..calculations.tkdict import FDFDict
from .base import SiestaBaseWorkChain

def accept_python_types(spec):
    '''
    Receives a ProcessSpec object and sets the serializer 
    for the inputs to `aiida.orm.nodes.data.base.to_aiida_type`.

    In this way, the user will be able to pass python objects, which
    will be automatically converted to aiida nodes.
    '''

    def serializer(obj):

        if isinstance(obj, list):
            return List(list=obj)
        else:
            return to_aiida_type(obj)

    spec.input = partial(spec.input, serializer=serializer)

class BaseIteratorWorkChain(WorkChain):
    '''
    General workflow that runs SIESTA simulations iteratively.

    The workflow itself can not be used. It relies in a method
    called `add_inputs` to modify the inputs at each iteration.
    See `ParameterIterator` and `AttributeIterator` for an example
    of this.
    '''

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Define the outline of the workflow, i.e. the order in which methods
        # should be executed
        spec.outline(
            cls.initialize,
            while_(cls.next_step)(
                cls.run_calc,
                cls.process_calc,
            ),
            cls.return_results
        )

        # Define all the inputs that this workchain expects

        # With this, we will automatically accept python datatypes, saving the user
        # the hassle of passing aiida types
        accept_python_types(spec)

        # Inputs related to the variable parameter
        spec.input(
            "init_value",
            valid_type=(Int, Float),
            required=False,
            help="The inital value of the parameter."
        )
        spec.input(
            "step",
            valid_type=(Int, Float),
            required=False,
            help="The step to apply for each iteration."
        )
        spec.input(
            'steps',
            valid_type=Int,
            required=False,
            help="The maximum number of steps to take"
        )
        spec.input(
            "values_list",
            valid_type=List,
            required=False,
            help='''The list of values to try. Use this if `init_value` and `step` are not
            suitable for you. E.g. the values are strings or the steps are not regular.'''
        )

        # Inputs related to the parameter to converge
        spec.input(
            "target",
            valid_type=Str,
            required=False,
            default=lambda: Str('E_KS'),
            help="The parameter that you want to track."
        )
        spec.expose_inputs(SiestaBaseWorkChain, exclude=('metadata',))
        spec.inputs._ports['pseudos'].dynamic = True
        
        # Define what are the outputs that this workchain will return
        spec.output('attempted_values')

    def initialize(self):

        self.ctx.variable_values = []
        self.ctx.target_values = []

        self.ctx.values_iterable = self._get_iterable()

    def _get_iterable(self):

        # The values to try will be provided by an iterator
        if 'values_list' in self.inputs:
            iterable = iter(self.inputs.values_list.get_list())
        elif hasattr(self, '_default_values_list') and 'init_value' not in self.inputs:
            iterable = iter(self._default_values_list)
        else:
            init_val = self.inputs.init_value.value if 'init_value' in self.inputs else self._default_init
            step = self.inputs.step.value if 'step' in self.inputs else self._default_step

            iterable = itertools.count(init_val, step) 

        return iterable    

    def next_step(self):
        '''
        Puts the next step in context.

        Returns
        ---------
        boolean
            Whether it should move to the next step or not (see this workchain's
            outline in the `define` method)
        '''

        # Give the oportunity to abort next cycle
        if not self.should_proceed:
            return False

        # Check if we already did the number of steps requested
        # Remember that we've not yet appended the next value, that's
        # why we check >= steps and not > steps.
        if 'steps' in self.inputs and len(self.ctx.variable_values) >= self.inputs.steps.value:
            return False

        # Otherwise, try to get a new value
        try:
            self.ctx.variable_values.append(self._next_val())

            self.report(f'Next value: {self.ctx.variable_values[-1]}')
        except StopIteration:
            # However, it's possible that there are no more values to try
            return False
        
        return True

    # The only logic BaseIteratorWorkChain knows is to run until there are no more
    # values to iterate or we've reached the number of steps
    # This method may be overwritten by child classes (see ConvergenceWorkChain)
    should_proceed = True

    def _next_val(self):
        '''
        Gets the next value to try.

        NOTE: Calling this method irreversibly 'outdates' the current value. Therefore
        it makes no sense to call it outside the `next_step` method.
        '''
        return next(self.ctx.values_iterable)

    @property
    def current_val(self):
        return self.ctx.variable_values[-1]

    def run_calc(self):
        '''
        Runs the calculation with the current value of the variable parameter.
        '''

        # Get the general inputs for the siesta run
        inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))
        # Let the
        self.add_inputs(inputs)
        
        # Run the SIESTA simulation and store the results
        calculation = self.submit(SiestaBaseWorkChain, **inputs)

        return ToContext(calculation=calculation)

    def process_calc(self):
        pass

    def return_results(self):
        '''
        Takes care of returning the results of the workchain to the user.
        '''
        pass

class ParameterIterator(BaseIteratorWorkChain):

    _units = None

    @classmethod
    def define(self, spec):
        super().define(spec)

        # If the user is using directly this class, give the option to choose
        # the parameter they want to converge
        if not hasattr(self, '_parameter'):
            spec.input(
                'parameter',
                valid_type=Str,
                required=True,
                help='The parameter that you want to vary to find convergence'
            )

            spec.input('units', valid_type=Str, required=False, help='The units of the parameter')
    
    def initialize(self):
        super().initialize()

        if hasattr(self, '_parameter'):
            self.ctx.parameter = self._parameter
        else:
            self.ctx.parameter = self.inputs.parameter.value

    def add_inputs(self, inputs):
        '''
        Adds the parameter with the appropiate value to the inputs that go
        to the SIESTA calculation.
        '''

        # Maybe we want to update some parameters from a dict different
        # than inputs.parameters(i.e. the basis dict)
        attr = getattr(self, '_parameters_attribute', 'parameters')

        # Convert the inputs to an fdf dict to avoid duplicate keys
        parameters = getattr(inputs, attr, DataFactory('dict')())
        parameters = FDFDict(parameters.get_dict())

        val = self.current_val 

        units = getattr(self.inputs, 'units', self._units)
        if units is not None:
            val = f'{val} {getattr(units, "value", units)}'

        parameters[self.ctx.parameter] = val

        # And then just translate it again to a dict to use it for parameters
        new_parameters = DataFactory('dict')(dict={key: val for key, (val, _) in parameters._storage.items()})
        setattr(inputs, attr, new_parameters)

class BasisParameterIterator(ParameterIterator):
    _parameters_attribute = 'basis'

class AttributeIterator(BaseIteratorWorkChain):

    @classmethod
    def define(self, spec):
        super().define(spec)

        # If the user is using directly this class, give the option to choose
        # the parameter they want to converge
        if not hasattr(self, '_attribute'):
            # We are going to add here inputs specific to converging a parameter 
            spec.input(
                'attribute',
                valid_type=Str,
                required=True,
                help='The attribute that you want to vary to find convergence'
            )

    def initialize(self):
        super().initialize()

        if hasattr(self, '_attribute'):
            self.ctx.attribute = self._attribute
        else:
            self.ctx.attribute = self.inputs.attribute.value

    def add_inputs(self, inputs):
        '''
        Adds the attribute with the appropiate value to the inputs that go
        to the SIESTA calculation.
        '''

        setattr(inputs, self.ctx.attribute, self.current_val)

class BasisSizeIterator(BasisParameterIterator):

    _parameter = 'pao-basissize'
    _default_values_list = ['SZ', 'SZP', 'DZ', 'DZP', 'TZ', 'TZP']

class MeshCutoffIterator(ParameterIterator):

    _parameter = 'meshcutoff'
    _units = 'Ry'
    _default_init = 100
    _default_step = 100

class KpointComponentIterator(AttributeIterator):

    _attribute = 'kpoints'

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Inputs related to the variable parameter
        spec.input(
            "component",
            valid_type=Int,
            required=False,
            help="The k component to converge. One of {0,1,2}"
        )
    
    @property
    def current_val(self):

        k_val = super().current_val

        # If there is already some kpoints, then we are going to use them
        if 'kpoints' in self.inputs:
            mesh, offset = self.inputs.kpoints.get_kpoints_mesh()
        else:
            # Otherwise we will just create it
            mesh = [1,1,1]
            offset = [0,0,0]

        mesh[self.inputs.component.value] = int(k_val)

        k_mesh = KpointsData()
        k_mesh.set_kpoints_mesh(mesh, offset)

        return k_mesh

class EnergyShiftIterator(BasisParameterIterator):

    _parameter = 'pao-energyshift'
    _units = 'Ry'


basis_params = FDFDict(
    paobasissize={
        'defaults': {
            'values_list': ['SZ', 'SZP', 'DZ', 'DZP', 'TZ', 'TZP']
        }
    },
    paoenergyshift={
        'defaults': {
            'units': 'Ry'
        }
    }
)

params = FDFDict(
    meshcutoff={
        'defaults': {
            'units': 'Ry',
            'init_value': 100,
            'step': 100
        }
    }
)

attributes = {}

preset_iters = {
    ParameterIterator: params,
    BasisParameterIterator: basis_params,
    AttributeIterator: attributes
}

def get_iterator_and_defaults(iterate_over):

    for IteratorClass, presets in preset_iters.items():

        # This try-except loop is a workaround because FDF dict
        # does not have regular dict behavior.
        # https://github.com/albgar/aiida_siesta_plugin/issues/47
        try:
            preset = presets[iterate_over]
            if preset is None:
                raise KeyError

            defaults = preset.get("defaults", {})
            break

        except KeyError:
            continue
    else:
        # If we haven't found it a default for it, we will try to
        # assume it
        if FDFDict.translate_key(iterate_over)[:3] == 'pao':
            IteratorClass = BasisParameterIterator
        else:
            IteratorClass = ParameterIterator

        defaults = {}

    # Add the setting that will indicate what are we iterating over
    if IteratorClass in (ParameterIterator, BasisParameterIterator):
        defaults['parameter'] = iterate_over
    elif IteratorClass is AttributeIterator:
        defaults['attribute'] = iterate_over
        
    return IteratorClass, defaults

def iterate(over, call_method='run', **kwargs):
    from aiida.engine import run, submit

    Iterator, defaults = get_iterator_and_defaults(over)

    execute = run if call_method is 'run' else submit

    return execute(Iterator, *args, **{**defaults, **kwargs})
