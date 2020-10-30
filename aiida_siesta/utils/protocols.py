import os
import yaml
from aiida.orm import Group
from aiida.common import exceptions


class ProtocolManager:
    """
    This class is meant to become the central engine for the management of protocols.

    With the word "protocol" we mean a series of suggested inputs for AiiDA WorkChains that allow
    users to more easly authomatize their workflows. Even though this approach could be general, at
    the moment we only think about protocols in the context of DFT inputs (Siesta inputs in our case).
    The choice of the inputs of a DFT simulation should be carefully tested for any new system.
    Users must be aware of the limitations of using protocols, but, in theory, this platform could
    become the place where we collect the "knowledge" about Siesta inputs. We hope that, with time,
    more and more protocols might be added to cover in a robust way entire categories of materials.
    This is the very beginning of the development and, for the moment, this engine is only able
    to collect protocol information from external files (a default `protocols_registry.yaml`
    and a custum-protocols file that a user put in the 'AIIDA_SIESTA_PROTOCOLS' environment variable)
    and build all the "core inputs" of a siesta calculation starting from a structure (some
    element-specific additions are performed - we call them atom_heuristics) and the choice
    of having relaxation and spin options. In the future we will probably need some kind of protocol
    algebra for merging, overriding, etc. protocols. This "core inputs" are then used by <WorkChain>
    specific input generators to create a ready-to submit sets of inputs.
    This class is, in fact, parent class of the class InputsGenerator, the fundation block
    of a series of classes <WorkChain>InputsGenerator with the scope
    to facilitate the choice of inputs for the corresponding <WorkChain>.

    The management of the pseudos is, at the moment, very fragile. It imposes that the user
    loads a pseudo_family with the exact same name of the one hard-coded for the protocol.
    Some other methods are implemented with the scope to access information about protocols and
    the API to use the protocols.
    """

    def __init__(self):
        """
        Construct an instance of ProtocolManager, validating the class attribute _calc_types set by the sub class
        and the presence of correct sintax in the protocols files (custom protocols can be set by users).
        """
        self._initialize_protocols()

        #Here we chack that each protocols implement correct syntax and mandatory entries
        self._protocols_checks()

    def _initialize_protocols(self):
        filepath = os.path.join(os.path.dirname(__file__), 'protocols_registry.yaml')

        with open(filepath) as thefile:
            self._protocols = yaml.full_load(thefile)

        if 'AIIDA_SIESTA_PROTOCOLS' in os.environ:
            bisfilepath = os.environ['AIIDA_SIESTA_PROTOCOLS']
            try:
                with open(bisfilepath) as thefile:
                    custom_protocols = yaml.full_load(thefile)
            except (IsADirectoryError, FileNotFoundError):
                raise RuntimeError(
                    'The environment variable devoted to custom protocols (AIIDA_SIESTA_PROTOCOLS) is set '
                    'to a not existent file'
                )
            self._protocols = {**self._protocols, **custom_protocols}

        self._default_protocol = 'standard_psml'

    def _protocols_checks(self):
        """
        Here implemented all the checks on the correct structure of each protocol. It also checks
        that, for each protocol, the correct pseudo family already loaded in the database.
        """

        def raise_invalid(message):
            raise RuntimeError('invalid protocol registry `{}`: '.format(self.__class__.__name__) + message)

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
                        '`{}`. Value and units are required'.format(k)
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

            #Run through heuristics
            for kind in structure.kinds:
                need_to_apply = False
                try:
                    cust_param = atomic_heuristics[kind.symbol]["parameters"]
                    need_to_apply = True
                except KeyError:
                    pass
                if need_to_apply:
                    if 'mesh-cutoff' in cust_param:
                        try:
                            cust_meshcut = float(cust_param["mesh-cutoff"].split()[0])
                        except (ValueError, IndexError):
                            raise RuntimeError(
                                "Wrong `mesh-cutoff` value for heuristc "
                                "{0} of protocol {1}".format(kind.symbol, key)
                            )
                        if meshcut_glob is not None:
                            if cust_meshcut > float(meshcut_glob):
                                meshcut_glob = cust_meshcut
                        else:
                            meshcut_glob = cust_meshcut
                            try:
                                meshcut_units = cust_param["mesh-cutoff"].split()[1]
                            except (ValueError, IndexError):
                                raise RuntimeError(
                                    "Wrong `mesh-cutoff` units for heuristc "
                                    "{0} of protocol {1}".format(kind.symbol, key)
                                )

            if meshcut_glob is not None:
                parameters["mesh-cutoff"] = "{0} {1}".format(meshcut_glob, meshcut_units)

        return parameters

    def _add_spin_options(self, key, orig_param):
        """
        Add to the parameters dictionary some additional parameters, called
        only if a calculation with spin polarization is requested.
        """
        if "spin_additions" in self._protocols[key]:
            #Check if mesh-cutoff is defined in "relax_additions" and has correct sintax
            if "mesh-cutoff" in self._protocols[key]["spin_additions"]:
                try:
                    v = self._protocols[key]["spin_additions"]
                    float(v["mesh-cutoff"].split()[0])
                    str(v["mesh-cutoff"].split()[1])
                except (ValueError, IndexError):
                    raise RuntimeError(
                        'Wrong format of `mesh-cutoff` in `spin_additions` of protocol '
                        '`{}`. Value and units are required'.format(key)
                    )
            #Merege the two dictionary in case of keys that are present in both dictionaries, the value
            #of the second dictionary is stored!
            parameters = {**orig_param, **self._protocols[key]["spin_additions"]}
            #Reset back mesh-cutoff if it was bigger in the original set.
            if "mesh-cutoff" in orig_param:
                if float(orig_param["mesh-cutoff"].split()[0]) > float(parameters["mesh-cutoff"].split()[0]):
                    parameters["mesh-cutoff"] = orig_param["mesh-cutoff"]
        else:
            parameters = orig_param.copy()

        return parameters

    def _add_relaxation_options(self, key, orig_param):
        """
        Add to the parameters dictionary some additional parameters, called
        only if a relaxation is requested
        """
        if "relax_additions" in self._protocols[key]:
            #Check if mesh-cutoff is defined in "relax_additions" and has correct sintax
            if "mesh-cutoff" in self._protocols[key]["relax_additions"]:
                try:
                    v = self._protocols[key]["relax_additions"]
                    float(v["mesh-cutoff"].split()[0])
                    str(v["mesh-cutoff"].split()[1])
                except (ValueError, IndexError):
                    raise RuntimeError(
                        'Wrong format of `mesh-cutoff` in `relax_additions` of protocol '
                        '`{}`. Value and units are required'.format(key)
                    )
            #Merege the two dictionary in case of keys that are present in both dictionaries, the value
            #of the second dictionary is stored!
            parameters = {**orig_param, **self._protocols[key]["relax_additions"]}
            #Reset back mesh-cutoff if it was bigger in the original set.
            if "mesh-cutoff" in orig_param:
                if float(orig_param["mesh-cutoff"].split()[0]) > float(parameters["mesh-cutoff"].split()[0]):
                    parameters["mesh-cutoff"] = orig_param["mesh-cutoff"]
        else:
            parameters = orig_param.copy()

        return parameters

    def _get_basis(self, key, structure):  # noqa: MC0001  - is mccabe too complex funct -
        """
        Method to construct the `basis` input.
        Heuristics are applied, a dictionary with the basis is returned.
        """
        basis = self._protocols[key]["basis"].copy()

        if "atomic_heuristics" in self._protocols[key]:  #pylint: disable=too-many-nested-blocks
            atomic_heuristics = self._protocols[key]["atomic_heuristics"]

            pol_dict = {}
            size_dict = {}
            pao_block_dict = {}

            #Run through all the heuristics
            for kind in structure.kinds:
                need_to_apply = False
                try:
                    cust_basis = atomic_heuristics[kind.symbol]["basis"]
                    need_to_apply = True
                except KeyError:
                    pass
                if need_to_apply:
                    if 'split-tail-norm' in cust_basis:
                        basis["pao-split-tail-norm"] = True
                    if 'polarization' in cust_basis:
                        pol_dict[kind.name] = cust_basis['polarization']
                    if 'size' in cust_basis:
                        size_dict[kind.name] = cust_basis['size']
                    if 'pao-block' in cust_basis:
                        pao_block_dict[kind.name] = cust_basis['pao-block']
                        if kind.name != kind.symbol:
                            pao_block_dict[kind.name] = pao_block_dict[kind.name].replace(kind.symbol, kind.name)

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
            if pao_block_dict:
                card = '\n'
                for k, v in pao_block_dict.items():
                    card = card + '{0} \n'.format(v)
                card = card + '%endblock pao-basis'
                basis['%block pao-basis'] = card

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

    def _get_pseudos(self, key, structure):

        from aiida_siesta.data.common import get_pseudos_from_structure

        family = self._protocols[key]["pseudo_family"]
        pseudos = get_pseudos_from_structure(structure, family)
        return pseudos
