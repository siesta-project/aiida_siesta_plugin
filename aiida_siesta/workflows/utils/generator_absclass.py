from abc import ABCMeta, abstractmethod
from aiida_siesta.workflows.utils.protocols import ProtocolManager


class InputsGenerator(ProtocolManager, metaclass=ABCMeta):
    """
    This metaclass sets the structure for WorkChain specific inputs generators.
    The subclasses needs to define _calc_types (containing the schema for the
    computational resources to use), and the _workchain_class (the WorkChain the
    subclass will produce inputs for).
    Moreover the methods `get_inputs_dict` and `get_filled_builder` must be overidden.
    The `get_inputs_dict` will implement all the logic to construct all the inputs
    and return them in a dictionary. The `get_filled_builder`, in theory, could
    be general, but it calls `get_inputs_dict` and this method can have different
    signature when overridden.
    """

    _calc_types = None

    def __init__(self, workchain_class):
        """
        Construct an instance of ProtocolManager, validating the class attribute _calc_types set by the sub class
        and the presence of correct sintax in the protocols files (custom protocols can be set by users).
        """

        super().__init__()

        if self._calc_types is None:
            message = 'invalid inputs generator `{}`: does not define `_calc_types`'.format(self.__class__.__name__)
            raise RuntimeError(message)

        try:
            workchain_class.get_builder()
        except AttributeError:
            message = 'invalid inputs generator `{}`: the defined `_workchain_class` is not a valid process'.format(
                self.__class__.__name__
            )
            raise RuntimeError(message)

        self._workchain_class = workchain_class

    def how_to_pass_computation_options(self):
        message = (
            "Computational resources are passed to get_filled_builder with the argument "
            "`calc_engines`. It's a dictionary with the following structure:"
        )
        return message, self._calc_types  #self._calc_types  #.values()

    @abstractmethod
    def get_inputs_dict(self, structure, calc_engines, protocol, **kwargs):
        """
        Return a dictionary with all the inputs, according to a protocol. The dictionary
        must contain only keywords that the builder of the corresponding WorkChain accepts!
        I think we should allow to change signature of this method.
        """

    @abstractmethod
    def get_filled_builder(self, structure, calc_engines, protocol, **kwargs):
        """
        Here it's the place were one should call `get_inputs_dict` in order to
        obtain the dictionary of inputs and then call `_fill builder` to obtain the
        builder.
        I think we should allow to change signature of this method.
        """

    def _fill_builder(self, inp_dict):
        """
        Return a builder, prefilled. Needs _accepted_workchain to obtain the builder
        and in input `inp_dict`, the dictionary containing all the inputs
        """
        builder = self._workchain_class.get_builder()

        for k, v in inp_dict.items():
            if v:
                builder[k] = v

        return builder
