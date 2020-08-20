import itertools
from functools import partial
import numpy as np

from aiida.plugins import DataFactory
from aiida.engine import WorkChain, while_, ToContext
from aiida.orm import Str, List, Int, Node
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.orm.utils import load_node
from aiida.common import AttributeDict


class BaseIterator(WorkChain):
    '''
    General workflow that runs iteratively a given `_process_class`.

    This is an abstract class, but subclasses just need to define the attribute `_process_class`.

    The quantity to iterate over is defined in the input `iterate_over`. It is a dictionary ("key", "value")
    where "key" is the name of a parameter we want to iterate over (`str`) and "value" is a `list` with all
    the values to iterate over for the corresponding key. The `iterate_over` is a dictionary because it is
    possible to iterate over several keywords at the same time. The way the algorithm deals with these
    multiple iterations is decided by the `iterate_mode` input.
    Because aiida accepts in input only lists of json-serializable objects, we have a serializer
    that transforms the list of values into a list of pks of aiida objects storing those values.
    An example:
        struct1 = StructureData(ase=ase_struct_1)
        struct2 = StructureData(ase=ase_struct_2)
        iterate_over = {"structure" : [struct1,struct2]}
    will be internally serialized to iterate_over = {"structure" : [struct1.pk,struct2.pk]}.
    The serialization is managed in `_iterate_input_serializer` and `_values_list_serializer`, that can be
    overridden in case, but this will impose also the change of the `_next_val` method.

    By default this WorkChain allows iterations only directly over the inputs of `_process_class`, but it
    can be extended to support other parameters using the `_params_lookup` variable.
    This variable indicates how to parse the values when the user wants to iterate over a certain parameter.
    It should be a list of dictionaries where each dictionary should have the following keys:
        group_name: str, optional
            The name of the group of parameters. Currently not used, but it will probably be used
            to show help messages.
        input_key: str
            The input of the process where the parsed values will end up.
        parse_func: function
            The function that will be used to parse the values for the parameters. The first arguments
            of this function should be the value to parse (as an aiida node) and the full inputs
            AttributeDict that is going to be passed to the process class. Also it should accept the `input_key` and
            `parameter` keyword arguments, which provides the input_key where the parsed value will go
            and the name of the parameter that the function is parsing. Finally, it needs to accept all kwargs
            that you define in the `keys` key (see below).
            It should not modify the inputs AttributeDict in place, but return the parsed value.
            E.g.:
            def parse_myparameters(val, inputs, parameter, input_key, **kwargs):
                ...do your parsing
                return parsed_value
        condition: function, optional
            A function that receives the name of the parameter and decides whether it belongs to this group.
            The function should return either `True`, or `False`.
            It will only be used if the key is not explicitly defined in the `keys` key (see below).
            If not provided, it will always return `False`
        keys: dict
            A dictionary where each key is a parameter that can be accepted and each value is either
            a dict containing the kwargs that will be passed to `parse_func` or None (same as empty dict).
            Even when a parameter is not defined here, it will be accepted if it fulfills `condition`.
    The `_params_lookup` variable is only used in the `process_input_and_parse_func` method, so you can check
    the code there to understand how exactly is used (it is quite simple).
    IMPORTANT NOTE: The order in `_params_lookup` matters. The workchain will go group by group trying to
    match the input parameter. If it matches a certain group, it will settle with it and won't continue to
    check the following groups.

    The design of the class also makes this class extensible to do something more than just iterate.
    The methods cls.return_results, cls.analyze_process, cls._should_proceed can be overridden to support,
    for instance, a convergence check.

    Following, we display the expanded outline of the workchain in pseudocode to help
    understand how it works:
    -------------------------------------------------------------------------------
    cls.initialize:
      cls._parse_iterate_over: # sets up all the required internal variables
        for key in iterate_over.keys():
          cls.process_input_and_parse_func(key) # Here all the logic for supported iterations.
                                                # 1) The allowed "key" are defined.
                                                # 2) For each key a `_parsing_func` is stored. This function
                                                #   implements which input of _process_class needs to be modified
                                                #   and how in presence of this "key".
      cls._get_iterator
    while (cls.next_step): # cls._should_proceed and cls._store_next_val are called inside cls.next_step!
      cls.run_batch:
        for i in range(batch_size):
          if i!= 0:
            cls._store_next_val
          cls._run_process # The current value is under self.current_val and
                           # the keys under self.iteration_keys
      cls.analyze_batch:
        for process in batch:
          cls._analyze_process(process)
    cls.return_results
    --------------------------------------------------------------------------------
    '''

    # THE _process_class NEEDS TO BE PROVIDED IN CHILD CLASSES!
    _process_class = None

    def __init_subclass__(cls, *args, **kwargs):
        """
        Imposes some requirements for each class that will inherit from the present class
        """

        super().__init_subclass__(*args, **kwargs)

        process_class = getattr(cls, "_process_class", None)

        from aiida.engine import Process  # pylint: disable=cyclic-import

        if process_class is not None:
            try:
                if not issubclass(process_class, Process):
                    raise ValueError('no valid Process class defined for `_process_class` attribute')
            except TypeError:
                raise ValueError('no valid Process class defined for `_process_class` attribute')

    def __init__(self, *args, **kwargs):
        """
        Construct the instance. Just needed to implement a check on the passed _process_class
        """

        process_class = getattr(self, "_process_class", None)

        if process_class is None:
            raise ValueError('Trying to instanciate an Iterator with `_process_class` attribute not implemented')

        super().__init__(*args, **kwargs)

    # The _params_lookup variable contains the indications on allowed parameters to
    # iterate over. It is optional in child classes. If not set, only iteration
    # over the inputs of _process_class are allowed.
    _params_lookup = ()

    # The inputs of _process_class will always be exposed. This class variable is passed directly
    # to spec.expose_inputs and can be used, for instance, to exclude some inputs or defining namespaces.
    _expose_inputs_kwargs = {}  #{'exclude': ('metadata',), "namespace": 'process_class_inputs' }

    # Whether the created inputs should be reused by the next step instead of always
    # grabbing the initial user inputs (use case: sequential converger)
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
        if "namespace" in cls._expose_inputs_kwargs:
            cls._exposed_input_keys = spec._exposed_inputs[cls._expose_inputs_kwargs["namespace"]][cls._process_class]
            for input_key in cls._exposed_input_keys:
                spec.inputs._ports[cls._expose_inputs_kwargs["namespace"]][input_key].required = False
        else:
            cls._exposed_input_keys = spec._exposed_inputs[None][cls._process_class]
            for input_key in cls._exposed_input_keys:
                spec.inputs._ports[input_key].required = False

    def initialize(self):
        """
        Initializes the variables that are used through the workchain.
        The methods `_parse_iterate_over` and `_get_iterator` are called.
        The first method analyzes the input ports `iterate_over`, performs checks on its
        content and organizes into context the informations. The second method adds the info
        from the port `iterate_mode` and creates the iterator.
        """

        self.ctx.variable_values = []

        # Here we check the `iterate_over` input and prepare many variables that are used through the
        # workchain. They are:
        #  self.ctx.iteration_keys, list of keys of `iterate_over`
        #  self.ctx.iteration_vals, all the values for each key
        #  self.ctx._process_input_keys, dict containing a map between keys and inputs of _process_class
        #  self.ctx._parsing_funcs, dict containing a map between keys and functions to process them
        self._parse_iterate_over()

        # Here we create self.ctx.values_iterator, the actual iterator
        self.ctx.values_iterator = self._get_iterator()

    def _parse_iterate_over(self):
        """
        Sets up self.ctx.iteration_keys (list of keys of `iterate_over`, passed by user)
        and self.ctx.iteration_vals (list of lists. Each list contains the values
        for the corresponding key).
        For each element of self.ctx.iteration_keys, we run `process_input_and_parse_func`.
        It checks for forbidden keys and sets the management of allowed ones. A process_input
        (stored in self.ctx._process_input_key) and a parsing_func (in self.ctx._parsing_funcs)
        is associated to each key. See `process_input_and_parse_func` documentation for info.
        """

        iterate_over = self.inputs.iterate_over.get_dict()

        # Get the names of the parameters and the values
        self.ctx.iteration_keys = tuple(iterate_over.keys())
        self.ctx.iteration_vals = tuple(iterate_over[key] for key in self.ctx.iteration_keys)

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

        iterate_mode = self.inputs.iterate_mode.value

        # Define the iterator depending on the iterate_mode.
        if iterate_mode == 'zip':
            iterator = zip(*self.ctx.iteration_vals)
        elif iterate_mode == 'product':
            iterator = itertools.product(*self.ctx.iteration_vals)

        return iterator

    @classmethod
    def process_input_and_parse_func(cls, parameter):
        """
        Function that makes use of the _params_lookup attribute in order to understand
        accepted iteration parameters and their management.
        Scopes:
        1) To implement errors for forbidden parameters.
        2) For each allowed parameter, to define the process_input and the parsing_func. The process_input is the
           `_process_class` input that needs to be modified when the parameter is called. The parsing function is
           the function that will implement the modifications. Each parsing function must return a node
           of the type accepted by process_input.
           For instance if one wants to support iterations over `pao-energy-shift` of siesta, its process_input
           will be `basis` and parsing_func will be a function that takes the value of pao-energy-shift
           and modifies the `basis` port accordingly.
        Both this tasks are done reading through _params_lookup.

        :param parameter: parameter (`str`) the user wants to iterate over, key of `iterate_over`
        :return: a `str` with the name of the process input that needs to be modified in presence of "key"
        :return: a python function that takes "key" and a corresponding value and returns an aiida dat node
                 of the type accepted by "key".
        """

        #In order to make use of cls._exposed_input_keys, we need to call first spec
        cls.spec()

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
            for parameters_group in cls._params_lookup:

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
