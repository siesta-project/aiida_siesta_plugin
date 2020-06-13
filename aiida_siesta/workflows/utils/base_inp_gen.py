from aiida_siesta.workflows.utils.protocols import ProtocolManager


class BaseWorkChainInputsGenerator(ProtocolManager):
    """
    This class has two main purposes:
    1) Provide a method (get_builder) that returns a builder for the WorkChain
       SiestaBaseWorkChain with pre-compiled inputs according to a protocol and
       some relaxation/bands options. This builder can be submitted to perform a Siesta
       scf, scf+bands, relaxation.
    2) Implement few methods that can be used to explore the relaxation types available,
       the bands options and protocols options (inherited by ProtocolRegistry)
    """

    _calc_types = {
        "siesta": {
            'code': 'Put here the code name, must be for plugin siesta.siesta',
            'options': 'Put here the computational options for running the relaxation, following the usual '
            'aiida schema'
        }
    }
    _relax_types = {
        "atoms_only":
        "the lattice shape and volume are fixed, only the atomic positions are relaxed",
        "variable_cell":
        "the lattice is relaxed together with the atomic coordinates. It allows "
        "to target hydro-static pressures or arbitrary stress tensors.",
        "constant_volume":
        "the cell volume is kept constant in a variable-cell relaxation: only "
        "the cell shape and the atomic coordinates are allowed to change.  Note that "
        "it does not make much sense to specify a target stress or pressure in this "
        "case, except for anisotropic (traceless) stresses"
    }
    _bands_path_generators = {
        "seekpath":
        "primitive cell always used. Always ok",
        "legacy":
        "To be used only if it is important to keep unchanged the cell shape"
        "of the input sturcture (no primitive cell). Not guaranteed to give correct paths",
    }

    def __init__(self, workchain_class):
        """
        Construct an instance of the inputs generator, validating the class attributes.
        """

        super().__init__()

        print(workchain_class.get_name())

        self._workchain_class = workchain_class

        def raise_invalid(message):
            raise RuntimeError('invalid protocol registry `{}`: '.format(self.__class__.__name__) + message)

        if self._calc_types is None:
            message = 'invalid inputs generator `{}`: does not define `_calc_types`'.format(self.__class__.__name__)
            raise RuntimeError(message)

        if self._relax_types is None:
            message = 'invalid inputs generator `{}`: does not define `_relax_types`'.format(self.__class__.__name__)
            raise RuntimeError(message)

        if self._bands_path_generators is None:
            message = 'invalid inputs generator `{}`: does not define `_bands_path_generators`'.format(
                self.__class__.__name__
            )
            raise RuntimeError(message)

    #Some methods to return info about options that the get_builder of this class
    #can obtain as optional input parameters.
    def get_relaxation_types(self):
        return list(self._relax_types.keys())

    def get_rel_type_info(self, key):
        try:
            return self._relax_types[key]
        except KeyError:
            raise ValueError("Wrong relaxation type: no relax_type with name {} implemented".format(key))

    def get_bands_path_generators(self):
        return list(self._bands_path_generators.keys())

    def get_bands_path_generator_info(self, key):
        try:
            return self._bands_path_generators[key]
        except KeyError:
            raise ValueError("Wrong path generator type: no bands_path_generator with name {} implemented".format(key))

    # noqa: MC0001  - is mccabe too complex funct -
    def get_inputs_dict(self, structure, calc_engines, protocol, bands_path_generator=None, relaxation_type=None):
        """
        Method return a dictionary with the inputs of a SiestaBaseWorkChain, obtained accordind
        to the protocol, structure, calc_engines passed by the user (optionally bands_path_generator
        and relaxation_typ as well)
        """

        from aiida.orm import (Dict, KpointsData)
        from aiida.orm import load_code
        from aiida.tools import get_explicit_kpoints_path

        #Checks
        if bands_path_generator:
            if bands_path_generator not in self.get_bands_path_generators():
                raise ValueError("No `bands_path_generator` with name {} implemented".format(bands_path_generator))
        if relaxation_type:
            if relaxation_type not in self.get_relaxation_types():
                raise ValueError("No `relaxation_type` with name {} implemented".format(relaxation_type))
        if not self.is_valid_protocol(protocol):
            import warnings
            defpro = self.get_default_protocol_name()
            warnings.warn("no protocol implemented with name `{0}`, using default `{1}`".format(protocol, defpro))
            protocol = defpro
        if 'siesta' not in calc_engines:
            raise ValueError("Wrong syntax in `calc_engines`. Check method `how_to_pass_computation_options`.")

        #protocol_dict = self.get_protocol(protocol)

        #Bandskpoints (might change structure)
        if bands_path_generator:
            if bands_path_generator == "legacy":
                legacy_kpath_parameters = {
                    'kpoint_distance': 0.05  # In units of b1, b2, b3 (Around 20 points per side...)
                }
                result = get_explicit_kpoints_path(structure, method='legacy', **legacy_kpath_parameters)
                ok_structure = structure
            else:
                seekpath_kpath_parameters = {'reference_distance': 0.01, 'symprec': 0.0001}
                result = get_explicit_kpoints_path(structure, **seekpath_kpath_parameters)
                ok_structure = result['primitive_structure']
            bandskpoints = KpointsData()
            bandskpoints = result['explicit_kpoints']
        else:
            ok_structure = structure
            bandskpoints = None

        #Parameters
        prot_param = self._get_param(protocol, ok_structure)
        #Add relaxation options in case requested
        if relaxation_type:
            parameters = self._add_relaxation_options(protocol, prot_param)
            parameters["md-type-of-run"] = "cg"
            parameters["md-num-cg-steps"] = 100
            if relaxation_type == "variable_cell":
                parameters["md-variable-cell"] = True
            if relaxation_type == "constant_volume":
                parameters["md-variable-cell"] = True
                parameters["md-constant-volume"] = True
        else:
            parameters = prot_param.copy()

        #Basis
        basis = self._get_basis(protocol, ok_structure)

        #Kpoints (might not be present, for molecules for instance)
        kpoints_mesh = self._get_kpoints(protocol, ok_structure)

        #Pseudo fam
        pseudo_fam = self._get_pseudo_fam(protocol)

        #Computational resources
        options = Dict(dict=calc_engines['siesta']["options"])
        code = load_code(calc_engines['siesta']["code"])

        inputs = {
            'structure': ok_structure,
            'parameters': Dict(dict=parameters),
            'code': code,
            'basis': Dict(dict=basis),
            'kpoints': kpoints_mesh,
            'pseudo_fam': pseudo_fam,
            'options': options,
            'bandskpoints': bandskpoints
        }

        return inputs

    def get_builder(self, structure, calc_engines, protocol, bands_path_generator=None, relaxation_type=None):

        #from aiida_siesta.workflows.base import SiestaBaseWorkChain

        inputs = self.get_inputs_dict(structure, calc_engines, protocol, bands_path_generator, relaxation_type)

        builder = self._workchain_class.get_builder()
        builder.structure = inputs["structure"]
        builder.basis = inputs["basis"]
        builder.parameters = inputs["parameters"]
        if inputs["kpoints"]:
            builder.kpoints = inputs["kpoints"]
        if inputs["bandskpoints"]:
            builder.bandskpoints = inputs["bandskpoints"]
        builder.pseudo_family = inputs["pseudo_fam"]
        builder.options = inputs["options"]
        builder.code = inputs["code"]

        return builder
