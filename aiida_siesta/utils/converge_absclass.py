import numpy as np

from aiida.engine import calcfunction
from aiida.orm import Float, Str, List, Bool, Int, load_node
from aiida.plugins import DataFactory

from .iterate_absclass import BaseIterator


@calcfunction
def generate_convergence_results(iteration_keys, used_values, target_values, converged, converged_index):
    '''
    Generates the final output of the convergence workflows.
    '''

    convergence_results = {
        'converged': Bool(converged),
    }

    if converged:
        #pylint: disable=unnecessary-comprehension
        converged_parameters = DataFactory('dict')(
            dict={key: val for key, val in zip(iteration_keys, used_values[converged_index.value])}
        )

        convergence_results['converged_parameters'] = converged_parameters
        convergence_results['converged_target_value'] = Float(target_values[converged_index.value])

    return convergence_results


class BasicConverger(BaseIterator):
    '''
    Plugin to add to an Iterator workchain to convert it to a convergence workflow.

    To use it, just build a class that inherits from this plugin and the iterator
    that you want to use.

    Examples
    -----------
    class MyConverger(BasicConverger, MyIterator):
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
            below_thresh = np.where(diffs < self.ctx.inputs.threshold.value)[0]

            converged = len(below_thresh) > 0  #pylint: disable=len-as-condition
            if converged:
                self.ctx.converged_index = below_thresh[0]

            self.report(
                'Convergence criterium: '
                '{0}; Last step diffs: {1}'.format(
                    self.ctx.inputs.threshold.value, diffs[-len(self.ctx.last_step_processes):]
                )
            )

        return converged

    def _should_proceed(self):
        return not self.converged

    def _analyze_process(self, process_node):
        """
        Takes the output from the process and stores the value of the target.
        """
        # Append the value of the target property for the last run
        results = process_node.outputs

        simulation_outputs = results.output_parameters.get_dict()

        self.ctx.target_values.append(simulation_outputs[self.ctx.inputs.target.value])

        super()._analyze_process(process_node)

    def return_results(self):
        '''
        Takes care of returning the results of the convergence to the user
        '''

        converged = Bool(self.converged)
        converged_index = Int(getattr(self.ctx, 'converged_index', -1))
        iteration_keys = List(list=list(self.ctx.iteration_keys))
        used_values = List(list=self.ctx.used_values)
        target_values = List(list=self.ctx.target_values)

        outputs = generate_convergence_results(iteration_keys, used_values, target_values, converged, converged_index)

        if converged:
            self.report(
                '\n\nConvergence has been reached! Converged parameters:'
                f'{outputs["converged_parameters"].get_dict()}\n'
            )
        else:
            self.report('\n\nWARNING: Workchain ended without finding convergence\n ')

        self.out_many(outputs)

        super().return_results()


class SequentialConverger(BaseIterator):
    '''
    Launches several convergence workchains sequentially. It is actally an iterator that
    iterate a _process_class which is a converger!!!!

    It is perfectly fine to define a class:
        class SiestaSeqConv(ProcessInputsIterator):
            _process_class = SiestaConverger
            _expose_inputs_kwargs = {"namespace": 'converger_inputs'}
    The presence of a namespace in _expose_inputs_kwargs is necessary to avoid conflict between the exposed
    iterate_over of _process_class  and the iterate_over of the SequentialConverger itself.
    However the class above will have a very convoluted way to submit calculations
        iterate_over = {"iterate_over" : [{"meshcutoff" : [...]},{"paoenergyshift":[....]}]
        inputs = { "converger_inputs" : {"inp1": ..., "inp2": ..., ... }  }
        submit(SiestaSeqConv, **inputs, iterate_over=iterate_over)
    and it will miss the part of incorporating the already converged "meshcutoff" into the convergence
    of "paoenergyshift".
    Here comes the need for a more complex class that facilitate the task of running sequential convergence.

    The iterate_over input is modified to accept a List where each item is a dict of accepted iterate_over
    keyword of _process_class (The Converger!!!). Imposes the modifications of the port and
    _iterate_input_serializer.

    THIS CLASS CAN NOT BE USED DIRECTLY, you need to subclass it and specify a cls._process_class,
    which should be the converger that you want to use.
    '''

    _expose_inputs_kwargs = {'exclude': ("iterate_over",), "namespace": 'converger_inputs'}
    _reuse_inputs = True

    def __init_subclass__(cls, *args, **kwargs):
        """
        Imposes some requirements for each class that will inherit from the present class
        """

        super().__init_subclass__(*args, **kwargs)

        if cls._process_class is not None:
            if not issubclass(cls._process_class, BasicConverger):
                raise ValueError(
                    'no valid Process class defined for `_process_class` attribute. It must be a Converger!'
                )

    @classmethod
    def define(cls, spec):
        super().define(spec)

        # In this case, iterate_over is going to be a list where each item is a dict as accepted
        # BaseIterator's regular iterate_over
        spec.inputs.ports['iterate_over'].valid_type = List

        # In principle we should not allow to provide a batch size and we should force a value of 1
        # But I don't know how to do it in a clean way. Also `iterate_mode` should be deeted from
        # the port, it is not used in this class.
        spec.inputs.ports['batch_size'].default = lambda: Int(1)

        spec.output(
            'converged_parameters',
            help="The values for the parameters that was enough to achieve convergence. "
            "If convergence is not achieved, it will be an empty dictionary",
        )
        spec.output('unconverged_parameters', required=False, help="The list of unconverged parameters.")

        spec.output('converged_target_value', required=False, help="The value of the target with convergence reached.")

    @classmethod
    def _iterate_input_serializer(cls, iterate_over):
        """
        Parses the "iterate_over" input of the workchain (a List of dicts) and applys serializations.
        For each dictionary in the list, each key is the name of a parameter we want to iterate
        over (str) and each value is a list with all the values to iterate over for that parameter.
        Therefore, for each element of the list, we just call the serializer of the of _process_class
        (That must be a converger!!!)

        :param iterate_over: aiida `List`. A list of dictionaries. See above.
        :return: an aiida Dict, the serialized input.
        """

        if isinstance(iterate_over, list):
            parsed_list = []
            for step in iterate_over:
                parsed_list.append(cls._process_class._iterate_input_serializer(step))
            return cls._values_list_serializer(parsed_list)

        return iterate_over

    def initialize(self):
        super().initialize()

        self.ctx.already_converged = {}

    def _parse_iterate_over(self):
        """
        The parameter to iterate over is `iterate_over` of _process_class (the converger!!!),
        the _process_input_keys corresponding to `iterate_over` is `iterate_over` itself and
        no parsing function is needed.
        """

        self.ctx.iteration_keys = ("iterate_over",)
        self.ctx._iteration_parsing = {"iterate_over": {}}

    def _get_iterator(self):
        """
        Iterator is now simply the list of values passed to `iterate_over`.
        """

        iterate_over = self.ctx.inputs.iterate_over.get_list()

        return zip(iterate_over)

    def _analyze_process(self, process_node):
        """
        Here we implement the part of incorporating the converged parameter in the
        sequential process.
        """

        if process_node.outputs.converged:

            # Get the parameters that have resulted in convergence
            new_converged = process_node.outputs.converged_parameters.get_dict()

            # Store them in what is going to be the final output
            self.ctx.already_converged.update(new_converged)

            # And now we are going to modify the inputs for the next step according
            # to what has converged, so that the convergence of the next parameter already incorporates
            # what we have found.
            # To modify the inputs we want to use _add_inputs of _process_class (the converger)
            # however this method makes use of self.ctx._iteration_parsing
            # of _process_class which are not the same for the present class.
            # Therefore we fake for a second this attribute.
            for key, val in new_converged.items():
                # Get the last inputs
                inputs = self.ctx.last_inputs
                proc_input, parse_func = self._process_class.process_input_and_parse_func(key)
                self.ctx._iteration_parsing[key] = {"input_key": proc_input, "parse_func": parse_func}
                # Modify "inputs" accordingly. Note that since we are using cls._reuse_inputs = True,
                # what we are modifying here will be used in the next iteration.
                self._process_class._add_inputs(self, key, val, inputs)

            self.ctx._iteration_parsing = {"iterate_over": {}}

            self.ctx.last_target_value = process_node.outputs.converged_target_value

        else:
            key = list(process_node.inputs.iterate_over.get_dict().keys())
            self.report(
                "\n WARINING: the next step of the sequential converger, if any, will be performed "
                f"with unchanged value for {key}, not with the last attempted value of the failed convergence."
            )

    def return_results(self):

        self.out('converged_parameters', DataFactory('dict')(dict=self.ctx.already_converged).store())

        #Create and return a list of parameters that did not converge. We compare
        #the iterate_over input and the self.ctx.already_converged. Both contains
        #untranslated keywords.
        iterate_over = self.ctx.inputs.iterate_over.get_list()
        # Get the list of all (possibly duplicated) parameters that have been iterated
        all_parameters = [key for pk_n in iterate_over for key in load_node(pk_n).get_dict()]
        # Find those that have not converged (while removing duplicates)
        unconv_list = list(set(all_parameters).difference(self.ctx.already_converged))
        if unconv_list:
            self.out('unconverged_parameters', DataFactory('list')(list=unconv_list).store())

        if "last_target_value" in self.ctx:
            self.out('converged_target_value', self.ctx.last_target_value)
