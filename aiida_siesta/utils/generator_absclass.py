from abc import ABCMeta, abstractmethod
from .protocols import ProtocolManager


class InputsGenerator(ProtocolManager, metaclass=ABCMeta):
    """
    This abstract class sets the structure for WorkChain specific inputs generators,
    meaning classes that are able to produce WorkChain inputs starting from
    a structure, a protocol and few more task specific options.

    Child classes need to define _calc_types (containing the schema for the
    computational resources to use) and the methods `get_inputs_dict` and `get_filled_builder`.
    The `get_inputs_dict` will implement all the logic to construct the inputs
    and return them in a dictionary. The `get_filled_builder` has the duty to create an already filled
    builder of the WorkChain. This last method, in theory, could be a general
    method that the child classes don't need to override. In fact it just call `get_inputs_dict`
    and places the dictionary entry in the empty builder.
    However the signature of `get_inputs_dict` can change case by case and the one of
    `get_filled_builder` must be change as well accordingly. If there was a possibility
    to impose automatically to `get_filled_builder` the same signature of `get_inputs_dict`,
    there woudn't be the need to define `get_filled_builder` as abstactmethod.

    Moreover, the workchain_class (the WorkChain the subclass will produce inputs for)
    is passed as a first argument during the instanciation of the class..
    This is not strictly necessary (the workchain_class could be requested as an attribute
    of child classes - exactly like _calc_types), however this strategy was adopted in order
    to avoid cyclic dependence and it can be advantageous if more then one workchain_class
    works with the exact same input generator.
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
        Return a builder, prefilled. Needs `_workchain_class` to obtain the builder
        and in input `inp_dict`, the dictionary containing all the inputs.
        """
        builder = self._workchain_class.get_builder()

        for k, v in inp_dict.items():
            if v is not None:
                builder[k] = v

        return builder
