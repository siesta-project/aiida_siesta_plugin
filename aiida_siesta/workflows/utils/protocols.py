import os
import yaml
from aiida.orm.groups import Group
from aiida.common import exceptions


class ProtocolRegistry:
    """
    This class is meant to become the central engine for the management of protocols.
    With the word "protocol" we mean a series of suggested inputs for AiiDA
    WorkChains that allow users to more easly authomatize their workflows.
    Even though this approach could be general, at the moment we only think
    about protocols in the context of DFT inputs (Siesta inputs in our case).
    The choice of the inputs of a DFT simulation should be carefully tested
    for any new system. Users must be aware of the limitations of using protocols,
    but, in theory, this platform could become the place where we collect
    the "know how" about Siesta inputs. We hope that, with time, more and more
    protocols might be added to cover in a robust way entire categories of materials.
    This is the very beginning of the development and, for the moment, only few
    methods and two very basic protocols are implemented.
    Moreover the methods are just functions useful to retrieve information about
    protocols, but in the future we will probably need some kind of protocol algebra
    for merging, overriding, etc. protocol values.
    The management of the pseudos is, in particular, very fragile. It imposes that the user
    loads a pseudo_family with the exact same name of the one hard-coded for the
    protocol.
    """

    filepath = os.path.join(os.path.dirname(__file__), 'protocols_registry.yaml')

    with open(filepath) as thefile:
        _protocols = yaml.full_load(thefile)

    _default_protocol = 'standard_delta'

    def __init__(self, *args, **kwargs):
        """Construct an instance of the protocol registry, validating the class attributes set by the sub class."""
        super().__init__(*args, **kwargs)

        def raise_invalid(message):
            raise RuntimeError('invalid protocol registry `{}`: '.format(self.__class__.__name__) + message)

        if not isinstance(self._protocols, dict):
            raise_invalid('protocols not collected in a dictionary')

        #Here we implement the checks on mandatory inputs we want for each protocol.
        #for the moment just a description
        for k, v in self._protocols.items():
            if not isinstance(self._protocols[k], dict):
                raise_invalid('protocol `{}` is not a dictionary'.format(k))

            if 'description' not in v:
                raise_invalid('protocol `{}` does not define the key `description`'.format(k))

            if 'pseudo_family' in v:
                famname = self._protocols[k]["pseudo_family"]
                try:
                    Group.get(label=famname)
                except exceptions.NotExistent:
                    raise_invalid(
                        'protocol `{}` requires `pseudo_family` with name {} '
                        'but no family with this name is loaded in the database'.format(k, famname)
                    )

        if self._default_protocol not in self._protocols:
            raise_invalid('default protocol `{}` is not a defined protocol'.format(self._default_protocol))

    def is_valid_protocol(self, name):
        return name in self._protocols

    def get_protocol_names(self):
        return list(self._protocols.keys())

    def get_default_protocol_name(self):
        return self._default_protocol

    def get_protocol_info(self, key):
        try:
            return self._protocols[key]["description"]
        except KeyError:
            raise ValueError("Wrong protocol: no protocol with name {} implemented".format(key))

    def get_protocol(self, key):
        try:
            return self._protocols[key]
        except KeyError:
            raise ValueError("Wrong protocol: no protocol with name {} implemented".format(key))
