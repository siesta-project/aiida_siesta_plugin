from abc import ABC
import itertools
from functools import partial
import numpy as np

from aiida.plugins import DataFactory
from aiida.engine import WorkChain, while_, ToContext
from aiida.orm import Str, List, KpointsData, Int, Node
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.orm.utils import load_node
from aiida.common import AttributeDict

from ...calculations.tkdict import FDFDict
from ..base import SiestaBaseWorkChain


class ProcessInputsIterator(WorkChain, ABC):
    '''
    General workflow that runs iteratively a given `_process_class`. This is an abstract class, but
    subclasses just need to define the attribute `_process_class`. By default this WorkChain allows iterations
    over all the inputs of `_process_class`, but it can be extended to support iterations over
    other parameters. The extention is done overriding the method `process_input_and_parse_func`.

    The quantity to iterate over is defined in the input `iterate_over`. It is a dictionary ("key", "value")
    where "key" is the name of a parameter we want to iterate over (`str`) and "value" is a `list` with all
    the values to iterate over for the corresponding key. The `iterate_over` is a dictionary because it is
    possible to iterate over several keywords at the same time. The way the algorithm deals with these
    multiple iterations is decided by the `iterate_mode` input.
    Because aiida accepts in input only lists of json-serializable objects, we have a serializer
    that transforms the list of values in a list of pks of aiida objects storing those values.
    An example:
        struct1 = StructureData(ase=ase_struct_1)
        struct2 = StructureData(ase=ase_struct_2)
        iterate_over = {"structure" : [struct1,struct2]}
    will be internally serialized to iterate_over = {"structure" : [struct1.pk,struct2.pk]}.
    The serialization is managed in `_iterate_input_serializer` and `_values_list_serializer`, that can be
    overridden in case, but this will impose also the change of the `_next_val` method.

    The design of the class also makes this class extensible to do something more than just iterate.
    The methods cls.return_results, cls.analyze_process, cls._should_proceed can be overridden to support,
    for instance, a convergence check.

    Following, we display the expanded outline of the workchain in pseudocode to help
    understand how it works and realise what can be overidden whithout breaking the code:
    -------------------------------------------------------------------------------
    cls.initialize:
        cls._get_iterator:
        for key in iterate_over.keys():
            cls.process_input_and_parse_func(key) # Here all the logic for supported iterations. Can be overridden.
                                                  # 1) The allowed "key" are defined.
                                                  # 2) For each key a `_parsing_func` is stored. This function
                                                  #   implements which input of _process_class needs to be modified
                                                  #   and how in presence of this "key".
    while (cls.next_step): # cls._should_proceed and cls._store_next_val are called inside cls.next_step!
        - cls.run_batch:
            for i in range(batch_size):
                if i!= 0:
                    cls._store_next_val
                cls._run_process # The current value is under self.current_val and
                                 # the keys under self.iteration_keys
        - cls.analyze_batch:
            for process in batch:
                cls.analyze_process(process)
    cls.return_results
    --------------------------------------------------------------------------------
    '''

    # THE PROCESS CLASS NEEDS TO BE PROVIDED IN CHILD CLASSES!
    _process_class = None

    def __init__(self, *args, **kwargs):
        """Construct the instance. Just needed to implement a check on the passed _process_class"""
        from aiida.engine import Process  # pylint: disable=cyclic-import

        try:
            if not issubclass(self._process_class, Process):
                raise ValueError('no valid Process class defined for `_process_class` attribute')
        except TypeError:
            raise ValueError('no valid Process class defined for `_process_class` attribute')

        super().__init__(*args, **kwargs)

    # The inputs of _process_class will always be exposed. This class variable is passed directly
    # to spec.expose_inputs and can be used, for instance, to exclude some inputs.
    _expose_inputs_kwargs = {'exclude': ('metadata',)}

    # Whether the created inputs should be reused by the next step
    # instead of always grabbing the initial user inputs (not sure use case)
    _reuse_inputs = False

    @classmethod
    def _iterate_input_serializer(cls, iterate_over):
        """
        Parses the "iterate_over" input of the workchain (a dictionary) and applys serializations.
        For each key-value of the dictionary, it parses the value (which is a list), using
        `cls._values_list_serializer`. See its documentation.

        :param iterate_over: `dict` or aiida `Dict`. A dictionary where each key is the name of a parameter
            we want to iterate over (`str`) and each value is a list with all the values to iterate over for
            that parameter.
        :return: an aiida Dict, the serialized input.
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
        Note that if you modify/overwrite this method you should take a look to the `_next_val` method.

        :param list_to_parse: a list with aiida data object
        :return: an aiida List containing the pk of each element of list_to_parse
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

        # Define the outline of the workflow, i.e. the order in which methods are executed.
        # See this class' documentation for an extended version of it
        spec.outline(cls.initialize,
                     while_(cls.next_step)(
                         cls.run_batch,
                         cls.analyze_batch,
                     ), cls.return_results)

        # Inputs that are general to all iterator workchains. They manage the iterator and batching.
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

        # We expose the inputs of the _process_class, in addition some more args
        # can be passed to the expose_inputs method (for instance inputs to exclude)
        spec.expose_inputs(cls._process_class, **cls._expose_inputs_kwargs)

        # Make the inputs that we have exposed not required (since if you iterate over,
        # them you might not pass them directly) and save them in a variable.
        cls._exposed_input_keys = spec._exposed_inputs[None][cls._process_class]

        for input_key in cls._exposed_input_keys:
            spec.inputs._ports[input_key].required = False

    def initialize(self):
        """
        Initializes the variables that are used through the workchain.
        The method `_get_iterator` is called. It analyzes the input ports
        `iterate_over` and `iterate_mode` and creates the iterator.
        The method process_input_and_parse_func is called. It defines the supported
        quantity to iterate over and their management.
        """

        self.ctx.variable_values = []

        # Here we create self.ctx.values_iterator, but also self.ctx.iteration_keys,
        # the list of keys passed by the user through the input `iterate_over`.
        self.ctx.values_iterator = self._get_iterator()

        # Here, for each key of `iterate_over`, we run `process_input_and_parse_func`.
        # It checks for forbidden keys and sets the management of allowed ones. A process_input
        # (stored in self.ctx._process_input_key) and a parsing_func (in self.ctx._parsing_funcs)
        # is associated to each key. See `process_input_and_parse_func` documentation for info.
        self.ctx._process_input_keys = {}
        self.ctx._parsing_funcs = {}
        for key in self.ctx.iteration_keys:
            proc_input, parsing_func = self.process_input_and_parse_func(key)

            self.ctx._process_input_keys[key] = proc_input
            self.ctx._parsing_funcs[key] = parsing_func

    def _get_iterator(self):
        """
        Builds the iterator that will generate values according to the inputs
        `iterate_over` and `iterate_mode`. Here also `port_name_and_parse_func` is called.
        """

        iterate_over = self.inputs.iterate_over.get_dict()
        iterate_mode = self.inputs.iterate_mode.value

        # Get the names of the parameters and the values
        self.ctx.iteration_keys = tuple(iterate_over.keys())
        self.ctx.iteration_vals = tuple(iterate_over[key] for key in self.ctx.iteration_keys)

        # Define the iterator depending on the iterate_mode.
        if iterate_mode == 'zip':
            iterator = zip(*self.ctx.iteration_vals)
        elif iterate_mode == 'product':
            iterator = itertools.product(*self.ctx.iteration_vals)

        return iterator

    @classmethod
    def process_input_and_parse_func(cls, key):  #pylint: disable=no-self-use
        """
        This method is an opportunity for the developers to define the accepted iteration parameters
        (keys of `iterate_over`) and their management. It is called in `initialize`.
        In the following, we will call "key" the parameter name (key of `iterate_over`).
        Scopes:
        1) To implement errors for forbidden keys.
        2) For each allowed key, to define the process_input and the parsing_func. The process_input is the
           `_process_class` input that needs to be modified when the key is called. The parsing function is
           the function that will implement the modifications. Each parsing function must return a node
           of the type accepted by process_input.
           For instance if one wants to support iterations over `meshcutoff` of siesta, its process_input
           will be `parameters` and parsing_func will be a function that takes the value of meshcutoff
           and modifies the parametrs port accordingly.
        By default, the only accepted keys are the name of the exposed inputs of `_process_class`, therefore
        the process_input is the key itself and there is no need for a parsing_function.

        :param key: parameter (`str`) the user wants to iterate over, key of `iterate_over`
        :return: a `str` with the name of the process input that needs to be modified in presence of "key"
        :return: a python function that takes "key" and a corresponding value and returns an aiida dat node
                 of the type accepted by "key".
        """

        if key not in cls._exposed_input_keys:
            raise ValueError("Iteration over {} is not supported".format(key))

        process_input = key
        parsing_func = None

        return process_input, parsing_func

    def next_step(self):
        """
        Puts the next step in context.
        Returns a boolean stating whether it should move to the next step or not (see this workchain's
        outline in the `define` method).
        """

        # Give the oportunity to abort next cycle
        if not self._should_proceed():
            return False

        # Otherwise, try to get a new value
        try:
            self._store_next_val()
        except StopIteration:
            # However, it's possible that there are no more values to try
            return False

        return True

    def _should_proceed(self):  #pylint: disable=no-self-use
        """
        The only logic a normal iterator knows is to run until there are no more
        values to iterate. However this method is an opportunity for developers to
        introduce a more complicated logic, for instance stop the iteration when a
        convergence is reached (see ConvergenceWorkChain).
        """
        return True

    def _next_val(self):
        '''
        Gets the next value to try. This method is called by `_store_next_val`.
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

    def _store_next_val(self):
        """
        Gets the next value from the iterator and appends it to the list of values.
        It also informs of what are the values that have been stored for use.
        """

        # Store the next value
        self.ctx.variable_values.append(self._next_val())

        # Inform about it
        info = '\n\t'.join([f'"{key}": {val}' for key, val in zip(self.ctx.iteration_keys, self.current_val)])
        self.report(f'Next values:{"{"}\n\t{info}\n{"}"}')

    @property
    def current_val(self):
        return self.ctx.variable_values[-1]

    def run_batch(self):
        '''
        Takes care of handling a batch of simulations. The numeber of processes for
        each batch is decided by the user through the input port `batch_size`.
        For each item in the batch, it calls `_run_process`.
        '''

        self.ctx.last_step_processes = []

        processes = {}
        batch_size = self.inputs.batch_size.value
        # Run as many processes as the "batch_size" input tells us to
        for i in range(batch_size):

            # If the batch size is bigger than 1, we need to retrieve more values
            if i != 0:
                try:
                    self._store_next_val()
                except StopIteration:
                    # But maybe there aren't enough values. In that case
                    # we will just run a smaller batch
                    break

            # Submit the process
            process_node = self._run_process()

            # Store the id of the process so that we can retrieve it later from context
            self.ctx.last_step_processes.append(process_node.uuid)

            # And then store the process (this will be passed to context at the end of the method)
            processes[process_node.uuid] = process_node

        self.report(f'Launched batch of {len(self.ctx.last_step_processes)}/{batch_size} processes')

        # Wait for the processes to finish
        return ToContext(**processes)

    def _run_process(self):
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
        # This is done like this so that _add_inputs is a simple method that
        # just takes a key and modifies the inputs accordingly. This should
        # be fine as most certainly each parameter can be set independently.
        # The self.current_val is the node containing the value!
        print(self.current_val)
        for key, val in zip(self.ctx.iteration_keys, self.current_val):
            self._add_inputs(key, val, inputs)

        self.ctx.last_inputs = inputs

        # Run the process and store the results
        process_node = self.submit(self._process_class, **inputs)

        return process_node

    def _add_inputs(self, key, val, inputs):
        '''
        Given a parameter (key) and the value (val, the node!), modify the inputs.
        We need:
        1) The `_process_class` input to modify (from `self.ctx._process_input_keys`)
        2) The parsed value, obtained applying the parsing function (from `self.ctx._parsing_funcs`)

        :param key: parameter (`str`) the user wants to iterate over, key of `iterate_over`
        :param val: an aiida data node containing the value associated to "key" at the present step
        :param inputs: AttributeDict containing all the `_process_class` inputs.
        '''

        # Get the name of the process input that will be modified
        attribute = self.ctx._process_input_keys.get(key, key)

        # Get the parsed value (it's the node), either the input value `val`
        # if no parsing_func is necessary or the value return by parsing_func.
        parsing_func = self.ctx._parsing_funcs.get(key, None)
        if parsing_func is not None:
            val = parsing_func(val, inputs)

        # Update the inputs AttributeDict
        setattr(inputs, attribute, val)

    def analyze_batch(self):
        '''
        Here, one could process the results of the processes that have been ran
        in the last batch.

        This default implementation just grabs each process in the order they have been
        launched and passes it to `analyze_process`. In this way, `analyze_process` does not need
        to handle anything related to the batch.
        '''

        for process_id in self.ctx.last_step_processes:

            self._analyze_process(self.ctx[process_id])

    def _analyze_process(self, process_node):
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


class GeneralIterator(ProcessInputsIterator):
    """
    This class proposes a general schema to extend the iteration over parameters other than the inputs
    of `_process_class`. It provides a systematic way to define groups of allowed parameters sharing
    the same process_input (the input of `_process_class` to be modified) and parsing_function
    (the function that returns a modifid process_input node according to the corresponding parameters).
    The core quantity of the schema is the attribute _params_lookup that is a list that must have
    the following  structure:
        _params_lookup = (
            ("group_name_1", {"input_keys": ..., "parse_func": ..., "condition": ..., "keys": ...})
            ("group_name_2", {"input_keys": ..., "parse_func": ..., "condition": ..., "keys": ...})
        )
    Explain here Pol, please

    ...
    It is an abstract class as its parent ProcessInputsIterator. The attributes _params_lookup
    and _process_class must be defined!

    """

    _params_lookup = ()

    def __init__(self, *args, **kwargs):
        """Construct the instance. Just needed to implement a check on the passed _params_lookup"""

        if not self._params_lookup:
            raise ValueError('The mandatory attribute `_params_lookup` is not defined')

        super().__init__(*args, **kwargs)

    @classmethod
    def process_input_and_parse_func(cls, parameter):  # pylint: disable=arguments-differ
        """
        Function that makes use of the _params_lookup attribute in order to understand
        accepted iteration parameters and their management.
        Scopes:
        1) To implement errors for forbidden keys.
        2) For each allowed key, to define the process_input and the parsing_func. The process_input is the
           `_process_class` input that needs to be modified when the key is called. The parsing function is
           the function that will implement the modifications. Each parsing function must return a node
           of the type accepted by process_input.
           For instance if one wants to support iterations over `meshcutoff` of siesta, its process_input
           will be `parameters` and parsing_func will be a function that takes the value of meshcutoff
           and modifies the parametrs port accordingly.
        Both this tasks are done reading through _params_lookup.

        :param parameter: parameter (`str`) the user wants to iterate over, key of `iterate_over`
        :return: a `str` with the name of the process input that needs to be modified in presence of "key"
        :return: a python function that takes "key" and a corresponding value and returns an aiida dat node
                 of the type accepted by "key".
        """

        # If the parameter is a key directly accepted by the process, that's it.
        # We will just interpret it as if the user wants to iterate over that input
        if parameter in cls._exposed_input_keys:
            input_key = parameter
            parse_func = None
            kwargs = {}

        # Otherwise, we will search for the parameter in the parameters look up list and
        # take the parse functions, input key and other arguments that are defined for it
        else:
            # Iterate through the parameters look up list until we find a group where
            # it belongs. Note that the parameter will stay in the first group that "accepts"
            # it, although there may be multiple groups where it fits. Therefore, order is important
            # in the parameters look up list
            for _, parameters_group in cls._params_lookup:

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


########################################################


class SiestaBaseWorkChainInputsIterator(ProcessInputsIterator):
    """
    Iterator for the SietaBaseWorkChain. No parameters other than the SietaBaseWorkChain inputs
    are allowed as keys of the input `iterate_over`.
    """
    _process_class = SiestaBaseWorkChain


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


class SieGenIterator(GeneralIterator):
    """
    Iterator for the SietaBaseWorkChain. The iterator is extended to iterate over any Siesta keyword.
    WARNING: if a keyword not recognized by Siesta is used in `iterate_over`, the iterator will not
    complain. It will just add the keyword to the parameters dict and run the calculation!
    """

    _process_class = SiestaBaseWorkChain
    _params_lookup = SIESTA_ITERATION_PARAMS
