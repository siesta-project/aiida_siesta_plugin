from aiida_siesta.workflows.utils.protocols import ProtocolRegistry


class SiestaBandsInputsGenerator(ProtocolRegistry):
    """
    This class is is similar to the SiestaRelaxInputsGenerator. It produces
    as well a builder for SiestaBaseWorkChain, but with correct inputs to perform
    a band structure analysis on a structure (without previous relaxation!!)

    """

    _calc_types = {
        "bands": {
            'code': 'Put here the code name, must be for plugin siesta.siesta',
            'options': 'Put here the computational options for running the relaxation, following the usual '
            'aiida schema'
        }
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

        if self._path_generators is None:
            message = 'invalid inputs generator `{}`: does not define `_path_generators`'.format(
                self.__class__.__name__
            )
            raise RuntimeError(message)

    def how_to_pass_computation_options(self):
        print(
            "The computational resources are passed to get_builder "
            "with the argument calc_engines. It is a dictionary with the "
            "following structure:"
        )
        return self._calc_types  #.values()

    def get_path_generators(self):
        return list(self._path_generators.keys())

    def get_path_generator_info(self, key):
        try:
            return self._path_generators[key]
        except KeyError:
            raise ValueError("Wrong path generator type: no path_generator with name {} implemented".format(key))

    def get_builder(self, structure, calc_engines, protocol, path_generator):

        from aiida_siesta.workflows.base import SiestaBaseWorkChain
        from aiida.orm import (Str, KpointsData, Dict)
        from aiida.orm import load_code
        from aiida.tools import get_explicit_kpoints_path

        #Checks
        #Checks
        if not self.is_valid_protocol(protocol):
            import warnings
            defpro = self.get_default_protocol_name()
            warnings.warn("no protocol implemented with name `{0}`, using default `{1}`".format(protocol, defpro))
            protocol = defpro
        if path_generator not in self.get_path_generators():
            raise ValueError(
                "Wrong path generator type: no path_generator with name {} implemented".format(path_generator)
            )
        if 'bands' not in calc_engines:
            raise ValueError("Wrong syntax in `calc_engines`. Check method `how_to_pass_computation_options`.")

        #Initialization
        protocol_dict = self.get_protocol(protocol)

        #Bandskpoints (might change structure)
        if path_generator == "legacy":
            legacy_kpath_parameters = {
                'kpoint_distance': 0.05  # In units of b1, b2, b3 (Around 20 points per side...)
            }
            result = get_explicit_kpoints_path(structure, method='legacy', **legacy_kpath_parameters)
            ok_structure = structure
        else:
            seekpath_kpath_parameters = {'reference_distance': 0.02, 'symprec': 0.0001}
            result = get_explicit_kpoints_path(structure, **seekpath_kpath_parameters)
            ok_structure = result['primitive_structure']
        bandskpoints = KpointsData()
        bandskpoints = result['explicit_kpoints']

        #Parameters
        parameters = self._get_param(protocol, ok_structure)

        #Basis
        basis = self._get_basis(protocol, ok_structure)

        #Kpoints
        if "kpoints" in protocol_dict:
            kpoints_mesh = KpointsData()
            kpoints_mesh.set_cell_from_structure(ok_structure)
            kp_dict = protocol_dict["kpoints"]
            if "offset" in kp_dict:
                kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"], offset=kp_dict["offset"])
            else:
                kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"])

        #Pseudo fam
        pseudo_fam = protocol_dict["pseudo_family"]

        #builder construction
        builder = SiestaBaseWorkChain.get_builder()
        builder.structure = ok_structure
        builder.basis = Dict(dict=basis)
        builder.parameters = Dict(dict=parameters)
        if "kpoints" in protocol_dict:
            builder.kpoints = kpoints_mesh
        builder.pseudo_family = Str(pseudo_fam)
        builder.options = Dict(dict=calc_engines['bands']["options"])
        builder.code = load_code(calc_engines['bands']["code"])
        builder.bandskpoints = bandskpoints

        return builder
