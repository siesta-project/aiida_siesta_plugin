import os
from abc import ABCMeta, abstractmethod
import yaml
from aiida.orm.groups import Group
from aiida.common import exceptions


class ProtocolManager(metaclass=ABCMeta):
    """
    Subclass should define _calc_types
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

    if 'AIIDA_SIESTA_PROTOCOLS' in os.environ:
        bisfilepath = os.environ['AIIDA_SIESTA_PROTOCOLS']
        try:
            with open(bisfilepath) as thefile:
                _custom_protocols = yaml.full_load(thefile)
        except (IsADirectoryError, FileNotFoundError):
            raise RuntimeError(
                'The environment variable devoted to custom protocols (AIIDA_SIESTA_PROTOCOLS) is set '
                'to a not existent file'
            )
        _protocols = {**_protocols, **_custom_protocols}

    _default_protocol = 'standard_delta'

    _calc_types = None

    def __init__(self):
        """
        Construct an instance of ProtocolManager, validating the class attribute _calc_types set by the sub class
        and the presence of correct sintax in the protocols files (can be set by user).
        """

        def raise_invalid(message):
            raise RuntimeError('invalid protocol registry `{}`: '.format(self.__class__.__name__) + message)

        #Here we implement the checks on mandatory inputs and structure we want for each protocol.
        if not isinstance(self._protocols, dict):
            raise_invalid('protocols not collected in a dictionary')

        for k, v in self._protocols.items():
            if not isinstance(self._protocols[k], dict):
                raise_invalid('protocol `{}` is not a dictionary'.format(k))

            if 'description' not in v:
                raise_invalid('protocol `{}` does not define the mandatory key `description`'.format(k))

            if 'parameters' not in v:
                raise_invalid('protocol `{}` does not define the mandatory key `parameters`'.format(k))
            if "mesh-cutoff" in v["parameters"]:
                try:
                    float(v["parameters"]["mesh-cutoff"].split()[0])
                    str(v["parameters"]["mesh-cutoff"].split()[1])
                except (ValueError, IndexError):
                    raise_invalid(
                        'Wrong format of `mesh-cutoff` in `parameters` of protocol '
                        '`{}`. Value and units required'.format(k)
                    )

            if 'basis' not in v:
                raise_invalid('protocol `{}` does not define the mandatory key `basis`'.format(k))

            if 'pseudo_family' not in v:
                raise_invalid('protocol `{}` does not define the mandatory key `pseudo_family`'.format(k))
            else:
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

        #here check that the subclass defines the attribute _calc_types
        if self._calc_types is None:
            message = 'invalid inputs generator `{}`: does not define `_calc_types`'.format(self.__class__.__name__)
            raise RuntimeError(message)

    #Some methods to return informations about the protocols
    #available and the _calc_types, describing the use of resources
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

    def how_to_pass_computation_options(self):
        print(
            "The computational resources are passed to get_builder with the "
            "argument `calc_engines`. It is a dictionary with the following structure:"
        )
        return self._calc_types  #.values()

    def _get_param(self, key, structure):
        """
        Method to construct the `parameters` input. Heuristics are applied, a dictionary
        with the parameters is returned.
        """
        parameters = self._protocols[key]["parameters"].copy()

        if "atomic_heuristics" in self._protocols[key]:  #pylint: disable=too-many-nested-blocks
            atomic_heuristics = self._protocols[key]["atomic_heuristics"]

            if 'mesh-cutoff' in parameters:
                meshcut_glob = parameters["mesh-cutoff"].split()[0]
                meshcut_units = parameters["mesh-cutoff"].split()[1]
            else:
                meshcut_glob = None

            for kind in structure.kinds:
                need_to_apply = False
                try:
                    cust_param = atomic_heuristics[kind.symbol]["parameters"]
                    need_to_apply = True
                except KeyError:
                    pass
                if need_to_apply:
                    if 'mesh-cutoff' in cust_param:
                        cust_meshcut = cust_param["mesh-cutoff"].split()[0]
                        if meshcut_glob:
                            if float(cust_meshcut) > float(meshcut_glob):
                                meshcut_glob = cust_meshcut
                        else:
                            meshcut_glob = cust_meshcut
                            meshcut_units = cust_param["mesh-cutoff"].split()[1]

            if meshcut_glob:
                parameters["mesh-cutoff"] = "{0} {1}".format(meshcut_glob, meshcut_units)

        return parameters

    def _add_relaxation_options(self, key, orig_param):
        if "relax" in self._protocols[key]:
            #in case of keys that are present in both dictionaries, the value
            #of the second dictionary is stored!
            parameters = {**orig_param, **self._protocols[key]["relax"]}
        else:
            parameters = orig_param.copy()

        return parameters

    def _get_basis(self, key, structure):  # noqa: MC0001  - is mccabe too complex funct -
        """
        Method to construct the `basis` input.
        Heuristics are applied, a dictionary with the basis is returned
        """
        basis = self._protocols[key]["basis"].copy()

        if "atomic_heuristics" in self._protocols[key]:  #pylint: disable=too-many-nested-blocks
            atomic_heuristics = self._protocols[key]["atomic_heuristics"]

            #Initializations
            #if 'pao-split-norm' in basis:
            #    split_norm_glob = basis["pao-split-norm"]
            #else:
            #    split_norm_glob = None
            #if 'pao-energy-shift' in basis:
            #    ene_shift_glob = basis["pao-energy-shift"].split()[0]
            #    en_shift_units = basis["pao-energy-shift"].split()[1]
            #else:
            #    ene_shift_glob = None
            pol_dict = {}
            size_dict = {}

            #Run through all the heuristics
            for kind in structure.kinds:
                need_to_apply = False
                try:
                    cust_basis = atomic_heuristics[kind.symbol]["basis"]
                    need_to_apply = True
                except KeyError:
                    pass
                if need_to_apply:
                    if 'split-norm' in cust_basis:
                        if cust_basis['split-norm'] == "tail":
                            basis["pao-split-tail-norm"] = True
                    #    else:
                    #        if kind.symbol == "H":
                    #            basis["pao-split-norm-h"] = cust_basis['split-norm']
                    #        else:
                    #            if split_norm_glob:
                    #                if float(cust_basis['split-norm']) > float(split_norm_glob):
                    #                    split_norm_glob = cust_basis['split-norm']
                    #            else:
                    #                split_norm_glob = cust_basis['split-norm']
                    if 'polarization' in cust_basis:
                        pol_dict[kind.name] = cust_basis['polarization']
                    #if 'energy-shift' in cust_basis:
                    #    cust_en_shift = cust_basis["energy-shift"].split()[0]
                    #cust_en_shif-units = cust_basis["energy-shift"].split()[1]
                    #    if ene_shift_glob:
                    #        if float(cust_en_shift) < float(ene_shift_glob):
                    #            ene_shift_glob = cust_en_shift
                    #    else:
                    #        ene_shift_glob = cust_en_shift
                    #        en_shift_units = cust_basis["energy-shift"].split()[1]
                    if 'size' in cust_basis:
                        size_dict[kind.name] = cust_basis['size']

            #Define the new basis dictionary
            #if split_norm_glob:
            #    basis["pao-split-norm"] = split_norm_glob
            #if ene_shift_glob:
            #    basis["pao-energy-shift"] = "{0} {1}".format(ene_shift_glob, en_shift_units)
            if pol_dict:
                card = '\n'
                for k, v in pol_dict.items():
                    card = card + '  {0}  {1} \n'.format(k, v)
                card = card + '%endblock paopolarizationscheme'
                basis['%block pao-polarization-scheme'] = card
            if size_dict:
                card = '\n'
                for k, v in size_dict.items():
                    card = card + '  {0}  {1} \n'.format(k, v)
                card = card + '%endblock paobasessizes'
                basis['%block pao-bases-sizes'] = card

        return basis

    def _get_kpoints(self, key, structure):
        from aiida.orm import KpointsData
        if "kpoints" in self._protocols[key]:
            kpoints_mesh = KpointsData()
            kpoints_mesh.set_cell_from_structure(structure)
            kp_dict = self._protocols[key]["kpoints"]
            if "offset" in kp_dict:
                kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"], offset=kp_dict["offset"])
            else:
                kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"])
        else:
            kpoints_mesh = None

        return kpoints_mesh

    def _get_pseudo_fam(self, key):
        from aiida.orm import Str
        return Str(self._protocols[key]["pseudo_family"])

    @abstractmethod
    def get_inputs_dict(self, structure, calc_engines, protocol):
        """
        Return a dictionary with all the inputs, according to protocol
        """

    @abstractmethod
    def get_builder(self, structure, calc_engines, protocol):
        """
        Return a builder, prefilled.
        """
