from abc import ABC, abstractmethod
import itertools
from functools import partial

from aiida.plugins import DataFactory
from aiida.engine import WorkChain, while_, ToContext, calcfunction
from aiida.orm import Float, Str, List, KpointsData, Int, Node, Dict
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.orm.utils import load_node
from aiida.common import AttributeDict

from ..calculations.tkdict import FDFDict
from .base import SiestaBaseWorkChain

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

    @classmethod
    def _iterate_input_serializer(cls, iterate_input):
        """
        Parses the input of an the iterator workchain.
        """

        if isinstance(iterate_input, dict):
            for key, val in iterate_input.items():
                iterate_input[key] = cls._values_list_serializer(val)
            
            iterate_input = DataFactory('dict')(dict=iterate_input)

        return iterate_input

    @staticmethod
    def _values_list_serializer(list_to_parse):
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

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Define the outline of the workflow, i.e. the order in which methods
        # should be executed
        spec.outline(cls.initialize,
                     while_(cls.next_step)(
                         cls.run_batch,
                         cls.analyze_batch,
                     ), cls.return_results)

        # Define all the inputs that this workchain expects

        # With this, we will automatically accept python datatypes, saving the user
        # the hassle of passing aiida types
        # accept_python_types(spec)

        # Inputs related to the variable parameter
        spec.input(
            "iterate_over",
            valid_type=DataFactory('dict'),
            serializer=cls._iterate_input_serializer, # The values list will always be a list of node pks.
            required=False,
            help='''A dictionary that indicates what to iterate over'''
        )
        spec.input(
            "iterate_mode", valid_type=Str, default=lambda: Str('zip'),
            help="Indicates the way the parameters should be iterated."
        )
        spec.input(
            "batch_size", valid_type=Int, default=lambda: Int(1),
            help="The number of simulations that should run at the same time"
        )

        # We expose the inputs of the workchain that is run at each iteration
        spec.expose_inputs(cls._process_class, exclude=cls._exclude_process_inputs, )
        spec.inputs._ports['pseudos'].dynamic = True

    def initialize(self):

        self.ctx.variable_values = []
        self.ctx.values_iterable = self._get_iterable()

    def _get_iterable(self):

        iterate_over = self.inputs.iterate_over.get_dict()
        iterate_mode = self.inputs.iterate_mode.value

        self.ctx.iteration_keys = tuple(iterate_over.keys())
        self.ctx.iteration_vals = tuple(iterate_over[key] for key in self.ctx.iteration_keys)
        self.ctx.iteration_keys = tuple(self._parse_key(key) for key in self.ctx.iteration_keys)

        if iterate_mode == 'zip':
            iterable = zip(*self.ctx.iteration_vals)
        elif iterate_mode == 'product':
            iterable = itertools.product(*self.ctx.iteration_vals)

        return iterable

    def _parse_key(self, key):
        return key

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

        # Otherwise, try to get a new value
        try:
            self.store_next_val()

        except StopIteration:
            # However, it's possible that there are no more values to try
            return False

        return True

    # The only logic BaseIteratorWorkChain knows is to run until there are no more
    # values to iterate or we've reached the number of steps
    # This method may be overwritten by child classes (see ConvergenceWorkChain)
    should_proceed = True

    def next_val(self):
        '''
        Gets the next value to try.

        NOTE: Calling this method irreversibly 'outdates' the current value. Therefore
        it makes no sense to call it outside the `next_step` method.
        '''

        next_pks = next(self.ctx.values_iterable)
        
        return tuple(load_node(next_pk) for next_pk in next_pks)
    
    def store_next_val(self):

        self.ctx.variable_values.append(self.next_val())

        info = '\n\t'.join([f'"{key}": {val}' for key, val in zip(self.ctx.iteration_keys, self.current_val)])

        self.report(f'Next values:{"{"}\n\t{info}\n{"}"}')
    
    @property
    def current_val(self):
        return self.ctx.variable_values[-1]

    def run_batch(self):
        '''
        Runs the calculation with the current value of the variable parameter.
        '''

        self.ctx.last_step_processes = []

        processes = {}
        batch_size = self.inputs.batch_size.value
        # Run as many processes as the "batch_size" input tells us to
        for i in range(batch_size):

            # If the batch size is bigger than 1, we need to retrieve more values
            if i != 0:
                try:
                    self.store_next_val()
                except StopIteration:
                    # But maybe there aren't enough values. In that case
                    # we will just run a smaller batch
                    break

            # Run the process and store the results
            process_node = self.run_process()

            self.ctx.last_step_processes.append(process_node.uuid)
            processes[process_node.uuid] = process_node

        self.report(f'Launched batch of {len(self.ctx.last_step_processes)}/{batch_size} processes')

        return ToContext(**processes)

    def run_process(self):

        # Get the general inputs for the siesta run
        inputs = AttributeDict(self.exposed_inputs(self._process_class))

        for key, val in zip(self.ctx.iteration_keys, self.current_val):
            self.add_inputs(key, val, inputs)

        # Run the process and store the results
        process_node = self.submit(self._process_class, **inputs)

        return process_node

    @abstractmethod
    def add_inputs(self, key, val, inputs):
        '''
        This method should be implemented in child classes.

        It takes all the inputs of the SIESTA calculation and it is expected to modify them
        in place according to the current iteration. See `ParameterIterator` and `AttributeIterator`
        for examples of this.
        '''

    def analyze_batch(self):
        '''
        Here, one could process the results of the simulation
        (which has been put into context)
        '''

        for process_id in self.ctx.last_step_processes:

            self.analyze_process(self.ctx[process_id])
    
    def analyze_process(self, process_node):
        """
        A child class has the oportunity to analyze a process here.
        """

    def return_results(self):
        '''
        Takes care of returning the results of the workchain to the user.
        '''


class AttributeIterator(BaseIteratorWorkChain):

    def _parse_val(self, key, val, inputs):
        return val

    def _attr_from_key(self, key):
        return key

    def add_inputs(self, key, val, inputs):
        '''
        Adds the attribute with the appropiate value to the inputs that go
        to the SIESTA calculation.
        '''

        attribute = self._attr_from_key(key)
        parsed_val = self._parse_val(key, val, inputs)

        setattr(inputs, attribute, parsed_val)


class SiestaIterator(AttributeIterator):

    _process_class = SiestaBaseWorkChain
    _exclude_inputs = ('metadata',)
    _attribute = None

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.input('units', valid_type=Str, required=False, help='The units of the parameter')

    def initialize(self):
        self.ctx._attributes = {}
        self.ctx._parsing_funcs = {}

        super().initialize()

    def _parse_key(self, key):
        
        attribute, parsing_func = self.attribute_and_parsing_func(key)

        self.ctx._attributes[key] = attribute
        self.ctx._parsing_funcs[key] = parsing_func

        return key
    
    def _attr_from_key(self, key):
        return self.ctx._attributes[key]

    def _parse_val(self, key, val, inputs):

        return self.ctx._parsing_funcs[key](val, inputs)

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
                attribute = param_info.get('attribute', attribute)
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

def set_up_kpoint_grid(val, inputs, attribute='kpoints', mode='distance'):

    old_kpoints = getattr(inputs, attribute, None)

    if old_kpoints is None:
        old_kpoints = KpointsData()

    try:
        mesh, offset = old_kpoints.get_kpoints_mesh()
    except (KeyError, AttributeError):
        mesh, offset = [1, 1, 1], [0, 0, 0]

    if not hasattr(old_kpoints, 'cell'):
        old_kpoints.set_cell_from_structure(inputs.structure)

    
    cell = old_kpoints.cell

    if mode == 'distance':
        new_kpoints = KpointsData()
        new_kpoints.set_cell(cell)
        new_kpoints.pbc = old_kpoints.pbc
        new_kpoints.set_kpoints_mesh_from_density(val.value, offset=offset)

    elif mode.endswith('component'):

        component = [ax for ax in range(3) if mode.startswith(str(ax))][0]

        mesh[component] = val

        new_kpoints = KpointsData()
        new_kpoints.set_kpoints_mesh(mesh, offset)

    return new_kpoints


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
        'kpoints': {},
        'kpoints_density':{
            'attribute': 'kpoints',
            'parse_func': set_up_kpoint_grid
        },
        **{
            f'kpoints_{ax}': {
                'attribute': 'kpoints',
                'parse_func': partial(set_up_kpoint_grid, mode=f'{ax}component')
            } for ax in range(3)
        }
    }
}
