from abc import ABC, abstractmethod
import itertools
from functools import partial

from aiida.plugins import DataFactory
from aiida.engine import WorkChain, while_, ToContext, calcfunction
from aiida.orm import Float, Str, List, KpointsData, Int, Node
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.orm.utils import load_node
from aiida.common import AttributeDict

from ..calculations.tkdict import FDFDict
from .base import SiestaBaseWorkChain

def normalize_list_to_pks(list_to_parse):
    '''
    Parses a list of objects to a list of node pks.

    This is done because aiida's List does not accept items with certain data structures
    (e.g. StructureData). In this way, we normalize the input to a list of pk, so that at
    each iteration we can access the value of the node.
    '''

    parsed_list = []
    for obj in list_to_parse:

        if not isinstance(obj, Node):
            obj = to_aiida_type(obj)

        if not obj.is_stored:
            obj.store()

        parsed_list.append(obj.pk)
    
    return List(list=parsed_list)

# def accept_python_types(spec):
#     '''
#     Receives a ProcessSpec object and sets the serializer
#     for the inputs to `aiida.orm.nodes.data.base.to_aiida_type`.

#     In this way, the user will be able to pass python objects, which
#     will be automatically converted to aiida nodes.
#     '''

#     def serializer(obj):

#         if isinstance(obj, list):
#             serialized = List(list=obj)
#         else:
#             serialized = to_aiida_type(obj)

#         return serialized

#     spec.input = partial(spec.input, serializer=serializer)


class BaseIteratorWorkChain(WorkChain, ABC):
    '''
    General workflow that runs simulations iteratively.

    The workchain itself can not be used.

    To use it, you need to define classes that inherit from it.
    These classes will have two main jobs:

        - Specify the calculation that needs to run iteratively.
    
    It relies in a method
    called `add_inputs` to modify the inputs at each iteration.
    See `ParameterIterator` and `AttributeIterator` for an example
    of this.
    '''

    # Here are some class variables that need to be defined
    @property
    @abstractmethod
    def _process_class(self):
        pass

    _exclude_process_inputs = ()
    _values_list_serializer = staticmethod(normalize_list_to_pks)

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Define the outline of the workflow, i.e. the order in which methods
        # should be executed
        spec.outline(cls.initialize,
                     while_(cls.next_step)(
                         cls.run_calc,
                         cls.process_calc,
                     ), cls.return_results)

        # Define all the inputs that this workchain expects

        # With this, we will automatically accept python datatypes, saving the user
        # the hassle of passing aiida types
        # accept_python_types(spec)

        # Inputs related to the variable parameter
        spec.input("init_value", valid_type=(Int, Float), required=False, help="The inital value of the parameter.")
        spec.input("step", valid_type=(Int, Float), required=False, help="The step to apply for each iteration.")
        spec.input('steps', valid_type=Int, required=False, help="The maximum number of steps to take")
        spec.input(
            "values_list",
            valid_type=List,
            serializer=cls._values_list_serializer, # The values list will always be a list of node pks.
            required=False,
            help='''The list of values to try. Use this if `init_value` and `step` are not
            suitable for you. E.g. the values are strings or the steps are not regular.'''
        )

        # We expose the inputs of the workchain that is run at each iteration
        spec.expose_inputs(cls._process_class, exclude=cls._exclude_process_inputs, )
        spec.inputs._ports['pseudos'].dynamic = True

    def initialize(self):

        self.ctx.variable_values = []
        self.ctx.values_iterable = self._get_iterable()

    def _get_iterable(self):

        # The values to try will be provided by an iterator
        if 'values_list' in self.inputs:
            iterable = iter(self.inputs.values_list.get_list())
        elif hasattr(self, '_default_values_list') and 'init_value' not in self.inputs:
            iterable = iter(self._default_values_list)
        else:
            try:
                init_val = self.inputs.init_value.value if 'init_value' in self.inputs else self._default_init
                step = self.inputs.step.value if 'step' in self.inputs else self._default_step
            except Exception:
                raise ValueError(
                    'You must specify either a "values_list" or an "init_value" and "step"'
                    ' so that we can build it for you'
                )

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

        next_pk = next(self.ctx.values_iterable)

        return load_node(next_pk)
    
    @property
    def current_val(self):
        return self.ctx.variable_values[-1]

    def run_calc(self):
        '''
        Runs the calculation with the current value of the variable parameter.
        '''

        # Get the general inputs for the siesta run
        inputs = AttributeDict(self.exposed_inputs(self._process_class))

        self.add_inputs(inputs)

        # Run the workchain and store the results
        workchain_futures = self.submit(self._process_class, **inputs)

        return ToContext(workchain_futures=workchain_futures)

    @abstractmethod
    def add_inputs(self, inputs):
        '''
        This method should be implemented in child classes.

        It takes all the inputs of the SIESTA calculation and it is expected to modify them
        in place according to the current iteration. See `ParameterIterator` and `AttributeIterator`
        for examples of this.
        '''

    def process_calc(self):
        '''
        Here, one could process the results of the simulation
        (which has been put into context)
        '''

    def return_results(self):
        '''
        Takes care of returning the results of the workchain to the user.
        '''

class AttributeIterator(BaseIteratorWorkChain):

    _attribute_required = True

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # If the user is using directly this class, give the option to choose
        # the parameter they want to converge
        if not hasattr(cls, '_attribute'):
            # We are going to add here inputs specific to converging a parameter
            spec.input(
                'attribute',
                valid_type=Str,
                required=cls._attribute_required,
                help='The attribute that you want to vary to find convergence'
            )

    def initialize(self, attribute=None):
        super().initialize()

        if attribute is None:
            if hasattr(self, '_attribute'):
                attribute = self._attribute
            else:
                attribute = self.inputs.attribute.value
        
        self.ctx.attribute = attribute

    def parse_val(self, val, inputs):
        return val

    def add_inputs(self, inputs):
        '''
        Adds the attribute with the appropiate value to the inputs that go
        to the SIESTA calculation.
        '''

        parsed_val = self.parse_val(self.current_val, inputs)

        setattr(inputs, self.ctx.attribute, parsed_val)

class SiestaIterator(AttributeIterator):

    _process_class = SiestaBaseWorkChain
    _exclude_inputs = ('metadata',)
    _attribute = None

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # If the user is using directly this class, give the option to choose
        # the parameter they want to converge
        if not hasattr(cls, '_parameter'):
            spec.input(
                'parameter',
                valid_type=Str,
                required=True,
                help='The parameter that you want to vary to find convergence'
            )

        spec.input('units', valid_type=Str, required=False, help='The units of the parameter')

    def initialize(self):

        if hasattr(self, '_parameter'):
            self.ctx.parameter = self._parameter
        else:
            self.ctx.parameter = self.inputs.parameter.value
        
        attribute, parsing_func = self.attribute_and_parsing_func(self.ctx.parameter)

        if parsing_func is not None:
            self.parse_val = parsing_func

        super().initialize(attribute=attribute)

    @staticmethod
    def attribute_and_parsing_func(parameter):
        '''
        Chooses which attribute to use for that convergence parameter
        '''

        if parameter in _INPUT_ATTRIBUTES['params']:
            attribute = parameter
            params_dict = _INPUT_ATTRIBUTES
        elif parameter.startswith('pao'):
            attribute = 'basis'
            params_dict = _BASIS_PARAMS
            parameter = FDFDict.translate_key(parameter)
        else:
            attribute = 'parameters'
            params_dict = _PARAMS
            parameter = FDFDict.translate_key(parameter)

        # Let's choose the parsing function
        if parameter in params_dict['params']:
            param_info = params_dict["params"].get(parameter, None)

            units = None
            try:
                units = param_info['defaults']['units']
            except:
                pass

            
            default_parse_func = params_dict.get('default_parse_func', None)
            if default_parse_func is not None:
                default_parse_func = partial(default_parse_func, parameter=parameter, units=units)

            if param_info is None:
                parsing_func = default_parse_func
            else:
                parsing_func = param_info.get('parse_func', default_parse_func)

        return attribute, parsing_func

def set_up_parameters_dict(val, inputs, parameter, attribute, units=None):

    val = val.value

    # Convert the inputs to an fdf dict to avoid duplicate keys
    parameters = getattr(inputs, attribute, DataFactory('dict')())
    parameters = FDFDict(parameters.get_dict())

    if units is not None:
        val = f'{val} {getattr(units, "value", units)}'

    parameters[parameter] = val

    # And then just translate it again to a dict to use it for parameters
    return DataFactory('dict')(dict={key: val for key, (val, _) in parameters._storage.items()})

class KpointComponentIterator(AttributeIterator):

    _attribute = 'kpoints'

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Inputs related to the variable parameter
        spec.input("component", valid_type=Int, required=False, help="The k component to converge. One of {0,1,2}")

    @property
    def current_val(self):

        k_val = super().current_val

        # If there is already some kpoints, then we are going to use them
        if 'kpoints' in self.inputs:
            mesh, offset = self.inputs.kpoints.get_kpoints_mesh()
        else:
            # Otherwise we will just create it
            mesh = [1, 1, 1]
            offset = [0, 0, 0]

        mesh[self.inputs.component.value] = int(k_val)

        k_mesh = KpointsData()
        k_mesh.set_kpoints_mesh(mesh, offset)

        return k_mesh

# Instead of defining classes for each param/attribute, another approach
# could be to have dicts defining "defaults" for each parameter attribute.
# For more complicated cases like kpoints we do need to define a new class though.

_BASIS_PARAMS = {
    'default_parse_func': partial(set_up_parameters_dict, attribute='basis'),
    'params': FDFDict(
        paobasissize={'defaults': {
            'values_list': ['SZ', 'SZP', 'DZ', 'DZP', 'TZ', 'TZP']
        }},
        paoenergyshift={'defaults': {
            'units': 'Ry'
        }}
    )
}

_PARAMS = {
    'default_parse_func': partial(set_up_parameters_dict, attribute='parameters'),
    'params': FDFDict(
        meshcutoff={'defaults': {'units': 'Ry', 'init_value': 100, 'step': 100}}
    )
}

# These are all the valid attributes for inputs (I believe)
_INPUT_ATTRIBUTES = {
    'default_parse_func': None,
    'params': {
        'code': {},
        'structure': {},
        'parameters': {},
        'pseudos': {},
        'basis': {},
        'settings': {},
        'parent_calc_folder': {},
    }
}


def iterate(over, call_method='run', **kwargs):
    '''
    Launches an aiida job to iterate through a parameter/basis parameter or attribute.

    Parameters
    -----------
    over: str
        the name of the parameter, basis parameter or attribute that you want to iterate
        over.
    call_method: {'run', 'submit'}, optional
        how the aiida job should be launched.
    **kwargs:
        All the extra keyword arguments that go directly into launching the job.
    '''
    from aiida.engine import run, submit

    execute = run if call_method == 'run' else submit

    return execute(SiestaIterator, **kwargs)
