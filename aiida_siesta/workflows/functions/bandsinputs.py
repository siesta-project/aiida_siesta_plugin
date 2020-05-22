from aiida_siesta.workflows.functions.protocols import ProtocolRegistry


class SiestaBandsInputsGenerator(ProtocolRegistry):
    """
    This class is is similar to the SiestaRelaxInputsGenerator. It produces
    as well a builder for SiestaBaseWorkChain, but with correct inputs to perform
    a band structure analysis on a structure (without previous relaxation!!)
    """

    _calc_types = {
        "bands": {
            'code_plugin': 'siesta.siesta',
            'description': 'These are calculations used for'
            'the main run of the code, running an scf and then the bands'
        }
    }
    _path_generators = {
        "seekpath":
        "primitive cell always used. Always ok",
        "legacy":
        "To be used only if it is important to keep unchanged the cell shape"
        "of the input sturcture (no primitive cell). Not guaranteed to give correct paths",
    }

    @classmethod
    def get_calc_types(cls):
        return cls._calc_types.keys()

    @classmethod
    def get_calc_type_schema(cls, key):
        try:
            return cls._calc_types[key]
        except KeyError:
            raise ValueError("Wrong calc_type: no calc_type {} implemented".format(key))

    @classmethod
    def get_path_generators(cls):
        return cls._path_generators.keys()

    @classmethod
    def get_builder(cls, structure, calc_engines, protocol, path_generator):

        from aiida_siesta.workflows.base import SiestaBaseWorkChain
        from aiida.orm import (Str, KpointsData, Dict)
        from aiida.orm import load_code
        from aiida.tools import get_explicit_kpoints_path

        if protocol not in cls.get_protocol_names():
            import warnings
            warnings.warn("no protocol implemented with name {}, using default standard".format(protocol))
            protocol = cls.get_default_protocol_name()

        protocol_dict = cls.get_protocol(protocol)

        #K points
        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(structure)
        kp_dict = protocol_dict["kpoints"]
        kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"], offset=kp_dict["offset"])

        #Parameters, including scf and relax options
        #scf
        parameters = protocol_dict["parameters"].copy()
        meshcutoff = 0
        min_meshcutoff = parameters["min_meshcut"]  # In Rydberg (!)
        del parameters["min_meshcut"]
        atomic_heuristics = protocol_dict["atomic_heuristics"]
        for kind in structure.get_kind_names():
            if atomic_heuristics[kind]:
                cutoff = atomic_heuristics[kind]['cutoff']
                meshcutoff = max(meshcutoff, cutoff)
        meshcutoff = max(min_meshcutoff, meshcutoff)
        parameters["meshcutoff"] = str(meshcutoff) + " Ry"
        #Walltime equal to scheduler, prevents calc to be killed by scheduler (it is implemented
        #in WorkChain as well, but this generator is more general
        parameters['max-walltime'] = calc_engines["bands"]["options"]["max_wallclock_seconds"]

        #Basis
        basis = protocol_dict["basis"]

        #Pseudo fam
        cls.is_pseudofamily_loaded(protocol)
        pseudo_fam = protocol_dict["pseudo_family"]

        #Bands
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

        builder = SiestaBaseWorkChain.get_builder()
        builder.structure = ok_structure
        builder.basis = Dict(dict=basis)
        builder.parameters = Dict(dict=parameters)
        builder.kpoints = kpoints_mesh
        builder.pseudo_family = Str(pseudo_fam)
        builder.options = Dict(dict=calc_engines["bands"]["options"])
        builder.code = load_code(calc_engines["bands"]["code"])
        builder.bandskpoints = bandskpoints

        return builder
