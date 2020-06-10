from aiida_siesta.workflows.utils.protocols import ProtocolRegistry


class BaseWorkChainInputsGenerator(ProtocolRegistry):
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
    _path_generators = {
        "seekpath":
        "primitive cell always used. Always ok",
        "legacy":
        "To be used only if it is important to keep unchanged the cell shape"
        "of the input sturcture (no primitive cell). Not guaranteed to give correct paths",
    }

    def __init__(self, *args, **kwargs):
        """Construct an instance of the inputs generator, validating the class attributes."""
        super().__init__(*args, **kwargs)

        def raise_invalid(message):
            raise RuntimeError('invalid protocol registry `{}`: '.format(self.__class__.__name__) + message)

        if self._calc_types is None:
            message = 'invalid inputs generator `{}`: does not define `_calc_types`'.format(self.__class__.__name__)
            raise RuntimeError(message)

        if self._relax_types is None:
            message = 'invalid inputs generator `{}`: does not define `_relax_types`'.format(self.__class__.__name__)
            raise RuntimeError(message)

        if self._path_generators is None:
            message = 'invalid inputs generator `{}`: does not define `_path_generators`'.format(
                self.__class__.__name__
            )
            raise RuntimeError(message)

    def how_to_pass_computation_options(self):
        print(
            "The computational resources are passed to get_builder with the "
            "argument `calc_engines`. It is a dictionary with the following structure:"
        )
        return self._calc_types  #.values()

    def get_relaxation_types(self):
        return list(self._relax_types.keys())

    def get_rel_type_info(self, key):
        try:
            return self._relax_types[key]
        except KeyError:
            raise ValueError("Wrong relaxation type: no relax_type with name {} implemented".format(key))

    def get_path_generators(self):
        return list(self._path_generators.keys())

    def get_path_generator_info(self, key):
        try:
            return self._path_generators[key]
        except KeyError:
            raise ValueError("Wrong path generator type: no path_generator with name {} implemented".format(key))

    # noqa: MC0001  - is mccabe too complex funct -
    def get_inputs(self, structure, calc_engines, protocol, path_generator=None, relaxation_type=None):  # noqa: MC0001

        from aiida.orm import (Str, KpointsData, Dict)
        from aiida.orm import load_code
        from aiida.tools import get_explicit_kpoints_path

        #Checks
        if path_generator:
            if path_generator not in self.get_path_generators():
                raise ValueError("No `path_generator` with name {} implemented".format(path_generator))
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

        protocol_dict = self.get_protocol(protocol)

        #Bandskpoints (might change structure)
        if path_generator:
            if path_generator == "legacy":
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
        parameters = self._get_param(protocol, ok_structure)
        #Add relaxation options in case requested
        if relaxation_type:
            parameters["md-type-of-run"] = "cg"
            parameters["md-num-cg-steps"] = 100
            if relaxation_type == "variable_cell":
                parameters["md-variable-cell"] = True
            if relaxation_type == "constant_volume":
                parameters["md-variable-cell"] = True
                parameters["md-constant-volume"] = True
            if "relax" in protocol_dict:
                parameters = {**parameters, **protocol_dict["relax"]}

        #Basis
        basis = self._get_basis(protocol, ok_structure)

        #Kpoints (might not be present, for molecules for instance)
        if "kpoints" in protocol_dict:
            kpoints_mesh = KpointsData()
            kpoints_mesh.set_cell_from_structure(ok_structure)
            kp_dict = protocol_dict["kpoints"]
            if "offset" in kp_dict:
                kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"], offset=kp_dict["offset"])
            else:
                kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"])
        else:
            kpoints_mesh = None

        #Pseudo fam
        pseudo_fam = Str(protocol_dict["pseudo_family"])

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

    def get_builder(self, structure, calc_engines, protocol, path_generator=None, relaxation_type=None):

        from aiida_siesta.workflows.base import SiestaBaseWorkChain

        inputs = self.get_inputs(structure, calc_engines, protocol, path_generator, relaxation_type)

        builder = SiestaBaseWorkChain.get_builder()
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
