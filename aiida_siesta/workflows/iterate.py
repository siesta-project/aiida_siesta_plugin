from abc import ABC, abstractmethod
import itertools
from functools import partial
import numpy as np

from aiida.plugins import DataFactory
from aiida.engine import WorkChain, while_, ToContext
from aiida.orm import Str, List, KpointsData, Int, Node
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.orm.utils import load_node
from aiida.common import AttributeDict

from ..calculations.tkdict import FDFDict
from .base import SiestaBaseWorkChain


class BaseIteratorWorkChain(WorkChain, ABC):
    '''
    General workflow that runs simulations iteratively.

    The workchain itself can not be used. To use it, you need to define classes
    that inherit from it. These classes will have one main job:

    Defining a `run_process` method that runs the process given the current
    value of the iterator. Note that this method does not need to take care of anything
    else! Following, we display the expanded outline of the workchain in pseudocode
    to help you understand how it works and realise what do you want to overwrite:

    -------------------------------------------------------------------------------
    cls.initialize:
        cls._get_iterator:
            keys = [self._parse_key(key) for key in keys]

    while (cls.next_step): # cls.should_proceed and cls.store_next_val are called inside cls.next_step!
        - cls.run_batch:
            for i in range(batch_size):
                if i!= 0:
                    cls.store_next_val
                cls.run_process # You have the current value under self.current_val and
                                # the keys under self.iteration_keys

        - cls.analyze_batch:
            for process in batch:
                cls.analyze_process(process)
    cls.return_results
    --------------------------------------------------------------------------------

    You have probably noted that you can analyze a process without worrying about the batch
    size by overwriting the `analyze_process` method.

    The main idea is that you probably don't need to subclass this class, because there is
    probably a subclass defined to serve a general purpose that your workflow fits into.
    See `InputIterator`, for example. The `BaseIteratorWorkchain` goal is just to provide
    a basis to handle the iterator and the batch size in a consistent way and facilitate
    further development.

    Finally, the class might seem unnecessarily complex due to the big number of methods
    defined in it. In reality, it is quite simple. The splitting of each step in separate
    methods is intentional to help the developer/user overwrite smalls bits of it without
    having to reimplement the whole step. Therefore: yes, you are invited to overwrite each
    little functionality to make it useful for your case :)

    You can even overwrite the input serializers. See `_iterate_input_serializer` and
    `_values_list_serializer`
    '''

    @classmethod
    def _iterate_input_serializer(cls, iterate_over):
        """
        Parses the "iterate_over" key of the workchain.

        For each key-value of the dictionary, it parses the value (which is a list),
        using `cls._values_list_serializer`. See its documentation.

        Parameters
        -----------
        iterate_over: dict or aiida Dict
            A dictionary where each key is the name of a parameter we want to iterate
            over (str) and each value is a list with all the values to iterate over for
            that parameter.

        Returns
        -----------
        aiida Dict
            the parsed input.
        """

        if isinstance(iterate_over, dict):
            for key, val in iterate_over.items():
                if not isinstance(val, (list, tuple, np.ndarray)):
                    raise ValueError(
                        f"We can not understand how to iterate over '{key}', "
                        f"you need to provide a list of values. You provided: {val}"
                    )
                iterate_over[key] = cls._values_list_serializer(val)

            iterate_over = DataFactory('dict')(dict=iterate_over)

        return iterate_over

    @staticmethod
    def _values_list_serializer(list_to_parse):
        '''
        Parses a list of objects to a list of node pks.

        This is done because aiida's List does not accept items with certain data structures
        (e.g. StructureData). In this way, we normalize the input to a list of pk, so that at
        each iteration we can access the value of the node.

        Note that if you modify/overwrite this method you should take a look to the `next_val` method.
        '''

        parsed_list = []
        # Let's iterate over all values so that we can parse them all
        for obj in list_to_parse:

            # If the object is a python type, convert it to a node
            if not isinstance(obj, Node):
                obj = to_aiida_type(obj)

            # If it has just been converted to a node, or it was an unstored node
            # store it so that it gets a pk.
            if not obj.is_stored:
                obj.store()

            # Now that we are sure the node has a pk, we append it to the list
            parsed_list.append(obj.pk)

        return List(list=parsed_list)

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Define the outline of the workflow, i.e. the order in which methods
        # should be executed. See this class' documentation for an extended
        # version of it
        spec.outline(cls.initialize,
                     while_(cls.next_step)(
                         cls.run_batch,
                         cls.analyze_batch,
                     ), cls.return_results)

        # Add inputs that are general to all iterator workchains. This is the stuff
        # that the BaseIteratorWorkchain manages for you: the iterator and batching.
        # More inputs can be defined in subclasses.
        spec.input(
            "iterate_over",
            valid_type=DataFactory('dict'),
            serializer=cls._iterate_input_serializer,
            help='''A dictionary where each key is the name of a parameter we want to iterate
            over (str) and each value is a list with all the values to iterate over for
            that parameter. Each value in the list can be either a node (unstored or stored)
            or a simple python object (str, float, int, bool).

            Note that each subclass might parse this keys and values differently, so you should
            know how they do it.
            '''
        )
        spec.input(
            "iterate_mode",
            valid_type=Str,
            default=lambda: Str('zip'),
            help='''Indicates the way the parameters should be iterated.
            Currently allowed values are:
                - 'zip': zips all the parameters together (all parameters should
                have the same number of values!)
                - 'product': performs a cartesian product of the parameters. That is,
                all possible combinations of parameters and values are explored.
             '''
        )
        spec.input(
            "batch_size",
            valid_type=Int,
            default=lambda: Int(1),
            help='''The maximum number of simulations that should run at the same time.
            You can set this to a very large number to make sure that all simulations run in
            one single batch if you want.'''
        )

    def initialize(self):
        """
        Initializes the variables that are used through the workchain.
        """

        self.ctx.variable_values = []
        self.ctx.values_iterator = self._get_iterator()

    def _get_iterator(self):
        """
        Builds the iterator that will generate values.
        """

        iterate_over = self.inputs.iterate_over.get_dict()
        iterate_mode = self.inputs.iterate_mode.value

        # Get the names of the parameters
        self.ctx.iteration_keys = tuple(iterate_over.keys())
        # Then the values
        self.ctx.iteration_vals = tuple(iterate_over[key] for key in self.ctx.iteration_keys)
        # And now that we've retrieved the values, parse the keys if necessary.
        # Note that by default _parse_key just returns the same key.
        self.ctx.iteration_keys = tuple(self._parse_key(key) for key in self.ctx.iteration_keys)

        # Define the iterator depending on the iterate_mode.
        if iterate_mode == 'zip':
            iterator = zip(*self.ctx.iteration_vals)
        elif iterate_mode == 'product':
            iterator = itertools.product(*self.ctx.iteration_vals)

        return iterator

    def _parse_key(self, key):  #pylint: disable=no-self-use
        """
        This method is an opportunity for the user to modify the iteration keys.

        It is called in `_get_iterator`
        """
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
        if not self.should_proceed():
            return False

        # Otherwise, try to get a new value
        try:
            self.store_next_val()
        except StopIteration:
            # However, it's possible that there are no more values to try
            return False

        return True

    # The only logic BaseIteratorWorkChain knows is to run until there are no more
    # values to iterate.
    # This method may be overwritten by child classes (see ConvergenceWorkChain)
    def should_proceed(self):  #pylint: disable=no-self-use
        return True

    def next_val(self):
        '''
        Gets the next value to try.

        This method is called by `store_next_val`.

        NOTE: Calling this method irreversibly 'outdates' the current value. Therefore
        it makes no sense to call it outside the `next_step` method.

        Since the input values are normalized to aiida pks, we load the corresponding
        nodes here. This method won't hold if the serializer is modified, so it should be
        changed accordingly.
        '''

        # Get the next values
        next_pks = next(self.ctx.values_iterator)

        # For each value, load the corresponding node
        return tuple(load_node(next_pk) for next_pk in next_pks)

    def store_next_val(self):
        """
        Gets the next value from the iterator and appends it to the list of values.

        It also informs of what are the values that have been stored for use.
        """

        # Store the next value
        self.ctx.variable_values.append(self.next_val())

        # Inform about it
        info = '\n\t'.join([f'"{key}": {val}' for key, val in zip(self.ctx.iteration_keys, self.current_val)])
        self.report(f'Next values:{"{"}\n\t{info}\n{"}"}')

    @property
    def current_val(self):
        return self.ctx.variable_values[-1]

    def run_batch(self):
        '''
        Takes care of handling a batch of simulations.

        For each item in the batch, it calls `run_process`.
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

            # Submit the process
            process_node = self.run_process()

            # Store the id of the process so that we can retrieve it later from context
            self.ctx.last_step_processes.append(process_node.uuid)

            # And then store the process (this will be passed to context at the end of the method)
            processes[process_node.uuid] = process_node

        self.report(f'Launched batch of {len(self.ctx.last_step_processes)}/{batch_size} processes')

        # Wait for the processes to finish
        return ToContext(**processes)

    @abstractmethod
    def run_process(self):
        """
        This method should be overwritten to implement running the process.
        """

    def analyze_batch(self):
        '''
        Here, one could process the results of the processes that have been ran
        in the last batch.

        This default implementation just grabs each process in the order they have been
        launched and passes it to `analyze_process`. In this way, `analyze_process` does not need
        to handle anything related to the batch.
        '''

        for process_id in self.ctx.last_step_processes:

            self.analyze_process(self.ctx[process_id])

    def analyze_process(self, process_node):
        """
        A child class has the oportunity to analyze a process here.

        Parameters
        -----------
        process_node:
            a process that has been ran in the batch.
        """

    def return_results(self):
        '''
        Takes care of postprocessing and returning the results of the workchain to the user.
        (if any outputs need to be returned)
        '''


class InputIterator(BaseIteratorWorkChain):
    """
    General workchain to iterate over a same process with different inputs.

    This class can not be used. However, you don't need to much!

    All you need to provide in a subclass to make use of this iterator is
    the process to run. The way you do this is by setting the `_process_class`
    variable. See the examples section . Once you have defined the process to run,
    this workchain exposes its inputs and runs the process iteritatively.

    The outline of this workchain is exactly the same as `BaseIteratorWorkchain`
    (see its docs). The only difference is that this class implements a `run_process`
    method that, in pseudocode, looks like this:

    cls.run_process:

        # Get the exposed inputs of the process to run
        inputs = get_exposed_inputs(MyProcess)

        # Modify the inputs
        for key, val in zip(iteration_keys, vals):
            cls.add_inputs(key, val, inputs):
                attribute = cls._attr_from_key(key)
                val = self._parse_val(key, val, inputs)
                setattr(inputs, attribute, parsed_val)

        # Submit the process with the modified inputs
        submit(MyProcess, **inputs)

    It also gives the possibility to accept iterating over parameters that are not
    directly an input of the process class.
    Therefore, if you want to add extra functionality to it, you probably should be
    good by overwriting the `_attr_from_key` or `_parse_val` methods, which by default
    do nothing. See `SiestaIterator` for an example of this.

    Examples
    ----------
    class MyIterator(InputIterator):
        _process_class = MyProcess

        # This class variable is passed directly to spec.expose_inputs
        _expose_inputs_kwargs = {'exclude': ('input_that_i_dont_want_to_expose', )}

    `MyIterator` will now run MyProcess iterating over all the input values.
    """

    # THE PROCESS CLASS NEEDS TO BE PROVIDED IN CHILD CLASSES!
    _process_class = None
    # The parser to use when a parameter is directly an accepted input of the process
    # class. See the 'input_key_and_parse_func' method to understand exactly how it is
    # used
    _process_parsing_func = None
    # The list where to look when a parameter that is not directly an input
    # of the process class is passed as parameter. See PARAMS at the bottom of
    # this file for an example of how it should look like.
    _params_lookup = ()
    # This class variable is passed directly to spec.expose_inputs
    _expose_inputs_kwargs = {'exclude': ('metadata',)}
    # Whether the created inputs should be reused by the next step
    # instead of always grabbing the initial user inputs
    _reuse_inputs = False

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # We expose the inputs of the workchain that is run at each iteration
        spec.expose_inputs(cls._process_class, **cls._expose_inputs_kwargs)

        # Get the list of inputs that the process class accepts.
        cls._process_input_keys = list(cls._process_class._spec.inputs._ports.keys())
        # And the ones that we have exposed.
        cls._exposed_input_keys = spec._exposed_inputs[None][cls._process_class]

        # Make the inputs that we have exposed not required (since if you iterate over,
        # them you might not pass them directly)
        for input_key in cls._exposed_input_keys:
            spec.inputs._ports[input_key].required = False

    def initialize(self):
        """
        Initialize the dicts to store input keys and parsing functions for
        each parameter.
        """
        self.ctx._input_keys = {}
        self.ctx._parsing_funcs = {}

        super().initialize()

    def run_process(self):
        '''
        Given a current value (self.current_val), runs the process.

        Before running, it sets up the inputs.
        '''

        # Get the exposed inputs for the process that we want to run.
        if self._reuse_inputs and hasattr(self.ctx, 'last_inputs'):
            inputs = self.ctx.last_inputs
        else:
            inputs = AttributeDict(
                self.exposed_inputs(self._process_class, namespace=self._expose_inputs_kwargs.get('namespace', None))
            )

        # For each key and value in this step, modify the inputs.
        # This is done like this so that add_inputs is a simple method that
        # just takes a key and modifies the inputs accordingly. This should
        # be fine as most certainly each parameter can be set independently.
        for key, val in zip(self.ctx.iteration_keys, self.current_val):
            self.add_inputs(key, val, inputs)

        self.ctx.last_inputs = inputs

        # Run the process and store the results
        process_node = self.submit(self._process_class, **inputs)

        return process_node

    def add_inputs(self, key, val, inputs):
        '''
        Given a parameter and the value, modify the inputs.

        We need to:
            - Get the input key (with `_attr_from_key`)
            - Get the parsed value (with `_parse_val`)

        Parameters
        -----------
        key: str
            the parameter that we need to set
        val: any
            the value for this parameter and this step.
        inputs: AttributeDict
            all the already existing inputs.
        '''

        # Get the input key
        attribute = self._attr_from_key(key)

        # Get the parsed value
        parsed_val = self._parse_val(val, key, inputs)

        # Update the inputs AttributeDict
        setattr(inputs, attribute, parsed_val)

    def _parse_key(self, key):
        """
        Uses the opportunity given by BaseIteratorWorkchain to store the input keys
        and parsing functions that we are going to use for each parameter through the workchain.

        Can be overwritten by child classes, but the default parameter management of this class
        is quite general and well organized, so consider giving it a try before implementing
        a different one.

        Parameters
        -----------
        key: str
            the name of the parameter, as passed in the "iterate_over" input.
        """

        input_key, parsing_func = self.input_key_and_parse_func(key)

        self.ctx._input_keys[key] = input_key
        self.ctx._parsing_funcs[key] = parsing_func

        return key

    def _parse_val(self, val, key, inputs):
        '''
        Parses the value using the parsing function that we have stored for this parameter.

        Can be overwritten by child classes, but the default parameter management of this class
        is quite general and well organized, so consider giving it a try before implementing
        a different one.

        Parameters
        -----------
        val: any
            the value to parse
        key: str
            the parameter to which this value belongs
        inputs: AttributeDict
            all the already existing inputs.
        '''

        parsing_func = self.ctx._parsing_funcs.get(key, None)

        if parsing_func is not None:
            val = parsing_func(val, inputs)

        return val

    def _attr_from_key(self, key):
        '''
        Returns the input key to which the value that corresponds to the parameter
        should be set.

        E.g.:
            if the parameter 'kpoints_whatever' means that you need to parse
            the value into a KpointsData() and set it to the 'kpoints' input,
            this should return 'kpoints'

        Can be overwritten by child classes, but the default parameter management of this class
        is quite general and well organized, so consider giving it a try before implementing
        a different one.

        Parameters
        ------------
        key: str
            the name of the parameter.
        '''
        return self.ctx._input_keys.get(key, key)

    def input_key_and_parse_func(self, parameter):
        '''
        Chooses which input key and parsing function to use for a parameter key.

        Parameters
        ----------
        parameter: str
            a key passed to the "iterate_over" input of the workchain.
        '''

        # Get the input keys accepted directly by the process
        process_input_keys = list(self._process_class._spec.inputs._ports.keys())

        # If the parameter is a key directly accepted by the process, that's it.
        # We will just interpret it as if the user wants to iterate over that input
        if parameter in process_input_keys:
            input_key = parameter
            parse_func = self._process_parsing_func
            kwargs = {}

        # Otherwise, we will search for the parameter in the parameters look up list
        # and take the parse functions, input key and other arguments that are defined for
        # it
        else:
            # Iterate through the parameters look up list until we find a group where
            # it belongs. Note that the parameter will stay in the first group that "accepts"
            # it, although there may be multiple groups where it fits. Therefore, order is important
            # in the parameters look up list
            for _, parameters_group in self._params_lookup:

                # Check if the parameter is in the group's explicit keys or if it matches a condition.
                in_keys = parameter in parameters_group['keys']
                if in_keys or parameters_group.get("condition", lambda p: False)(parameter):

                    # If it does, get the input key that this parameter is going to modify
                    input_key = parameters_group['input_key']
                    # The function that will parse the values before passing them to the process
                    parse_func = parameters_group['parse_func']
                    # And some extra arguments.
                    kwargs = parameters_group['keys'].get(parameter, None) or {}
                    break
            else:
                raise ValueError(f"We didn't find any known way to iterate over '{parameter}'")

        # Finally, make sure the parse function gets the input and parameter keys, as well as extra arguments.
        if parse_func is not None:
            parse_func = partial(parse_func, input_key=input_key, parameter=parameter, **kwargs)

        return input_key, parse_func


# The following are helper functions to parse input values in the SiestaIterator. See
# the global dict SIESTA_ITERATION_PARAMS to know which parameters make use of them.


def set_up_parameters_dict(val, inputs, parameter, input_key, defaults=None):
    """
    Parsing function that sets an fdf parameter.

    This is used by the basis parameters (input_key='basis') and the rest of
    fdf parameters (input_key='parameters')

    Parameters
    -----------
    val: str, float, int, bool
        the value to set for the parameter.
    inputs: AttributeDict
        all the current inputs, so that we can extract the curren FDFdict.
    parameter: str
        the name of the fdf flag. Case and hyphen insensitive.
    input_key: str
        the input of the fdf dict that should be modified.
    defaults: dict, optional
        a dictionary of defaults. The only key that matters right now is "units",
        which is the units to set for the value if the value is a number. This may change.
    """

    val = getattr(val, "value", val)

    # Get the current FDFdict for the corresponding input
    parameters = getattr(inputs, input_key, DataFactory('dict')())
    parameters = FDFDict(parameters.get_dict())

    # Set the units for the value if needed
    if isinstance(val, (int, float)) and defaults and "units" in defaults:
        val = f'{val} {defaults["units"]}'

    # Then set the value of the parameter in the FDF dict.
    parameters[parameter] = val

    # And then just translate it again to a dict to use it in the input
    return DataFactory('dict')(dict=parameters.get_dict())


def set_up_kpoint_grid(val, inputs, parameter, input_key='kpoints'):
    """
    Parsing function that sets a kpoint grid.

    This is used by the basis parameters (input_key='basis') and the rest of
    fdf parameters (input_key='parameters')

    Parameters
    -----------
    val: str, float, int, bool
        the value to set for the parameter.
    inputs: AttributeDict
        all the current inputs, so that we can extract the old kpoints.
    parameter: str, {'kpoints_density', 'kpoints_0', 'kpoints_1', 'kpoints_2'}
        used to understand how to parse the value.

        If 'kpoints_density': the value is interpreted as the maximum distance
        between two grid points along a reciprocal axis.

        Else, it is interpreted as the number of points for one of the components
        of the grid.
    input_key: str
        the input of the fdf dict that should be modified.
    """

    # If there is already a KpointsData() in inputs grab it to modify it
    old_kpoints = getattr(inputs, input_key, None)

    # Else define a new one
    if old_kpoints is None:
        old_kpoints = KpointsData()

    # Get the mesh and the offset
    try:
        mesh, offset = old_kpoints.get_kpoints_mesh()
    except (KeyError, AttributeError):
        mesh, offset = [1, 1, 1], [0, 0, 0]

    # Get also the cell
    if hasattr(old_kpoints, 'cell'):
        cell = old_kpoints.cell
        pbc = old_kpoints.pbc
    else:
        cell = inputs.structure.cell
        pbc = inputs.structure.pbc

    # And finally define the new KpointsData according to the required mode.

    # Change the density
    if parameter == 'kpoints_density':
        new_kpoints = KpointsData()
        new_kpoints.set_cell(cell)
        new_kpoints.pbc = pbc
        new_kpoints.set_kpoints_mesh_from_density(val.value, offset=offset)

    # Change an individual component
    else:

        component = int(parameter[-1])

        mesh[component] = val

        new_kpoints = KpointsData()
        new_kpoints.set_kpoints_mesh(mesh, offset)

    return new_kpoints


# This is the parameters' look up list for the siesta iterator, which enables iterating
# over extra parameters apart from the inputs of SiestaBaseWorkChain. This may be taken
# as an example to allow extra parameters in any input iterator that uses a different
# process.

SIESTA_ITERATION_PARAMS = (
    (
        "Basis parameters", {
            "input_key":
            "basis",
            "parse_func":
            set_up_parameters_dict,
            "condition":
            lambda parameter: FDFDict.translate_key(parameter).startswith("pao"),
            "keys":
            FDFDict({
                "paobasissize": {
                    'defaults': {
                        'values_list': ['SZ', 'SZP', 'DZ', 'DZP', 'TZ', 'TZP']
                    }
                },
                "paoenergyshift": {
                    'defaults': {
                        'units': 'Ry'
                    }
                }
            })
        }
    ),
    (
        "SCF Brillouin zone", {
            "input_key": "kpoints",
            "parse_func": set_up_kpoint_grid,
            "keys": {
                'kpoints_density': None,
                'kpoints_0': None,
                'kpoints_1': None,
                'kpoints_2': None
            }
        }
    ),
    (
        "FDF parameters", {
            "input_key": "parameters",
            "condition": lambda parameter: True,
            "parse_func": set_up_parameters_dict,
            "keys": FDFDict({"meshcutoff": {
                'defaults': {
                    'units': 'Ry',
                    'init_value': 100,
                    'step': 100
                }
            }})
        }
    ),
)


class SiestaIterator(InputIterator):
    """
    General iterator for ANY parameter in SIESTA.

    It runs `SiestaBaseWorkChain` iteratively setting the inputs
    appropiately at each iteration.

    With the help of the `SIESTA_ITERATION_PARAMS` dictionary, it decides
    what each parameter key in the "iterate_over" needs to do (i.e. which input
    to modify and how to parse the value into a valid input).

    Examples
    ------------

    See `examples/workflows/example_iterate.py`.
    """

    _process_class = SiestaBaseWorkChain
    _params_lookup = SIESTA_ITERATION_PARAMS
    _expose_inputs_kwargs = {'exclude': ('metadata',)}

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.inputs._ports['pseudos'].dynamic = True
        # spec.input('units', valid_type=Str, required=False, help='The units of the parameter')
