from aiida_siesta.workflows.utils.protocols import ProtocolRegistry


class SiestaRelaxInputsGenerator(ProtocolRegistry):
    """
    This class has two main purposes:
    1) Provide a method (get_builder) that returns a builder for the WorkChain
       SiestaBaseWorkChain with pre-compiled inputs according to a protocol and
       some relaxation options. This builder can be submitted to perform a Siesta
       relaxation.
    2) Implement few methods that can be used to explore the relaxation types available.
       (in the future they might be used to automatically generate a GUI that
       allow users to run a relaxation with Siesta after selecting options with
       few clicks).
    """

    _calc_types = {
        "relaxation": {
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

    def how_to_pass_computation_options(self):
        print(
            "The computational resources are passed to get_builder "
            "with the argument calc_engines. It is a dictionary with the "
            "following structure:"
        )
        return self._calc_types  #.values()

    def get_relaxation_types(self):
        return list(self._relax_types.keys())

    def get_rel_type_info(self, key):
        try:
            return self._relax_types[key]
        except KeyError:
            raise ValueError("Wrong relaxation type: no relax_type with name {} implemented".format(key))

    #pylint: disable=too-many-statements
    def get_builder(
        self, structure, calc_engines, protocol, relaxation_type, threshold_forces=None, threshold_stress=None
    ):

        from aiida_siesta.workflows.base import SiestaBaseWorkChain
        from aiida.orm import (Str, KpointsData, Dict)
        from aiida.orm import load_code

        #Checks
        if not self.is_valid_protocol(protocol):
            import warnings
            defpro = self.get_default_protocol_name()
            warnings.warn("no protocol implemented with name `{0}`, using default `{1}`".format(protocol, defpro))
            protocol = defpro
        if relaxation_type not in self.get_relaxation_types():
            raise ValueError("Wrong relaxation type: no relax_type with name {} implemented".format(relaxation_type))
        if 'relaxation' not in calc_engines:
            raise ValueError("Wrong syntax in `calc_engines`. Check method `how_to_pass_computation_options`.")

        protocol_dict = self.get_protocol(protocol)

        #Parameters
        parameters = self._get_param(protocol, structure)
        parameters["md-type-of-run"] = "cg"
        parameters["md-num-cg-steps"] = 100
        if relaxation_type == "variable_cell":
            parameters["md-variable-cell"] = True
        if relaxation_type == "constant_volume":
            parameters["md-variable-cell"] = True
            parameters["md-constant-volume"] = True
        if not threshold_forces:
            if "threshold-forces" in protocol_dict:
                threshold_forces = protocol_dict["threshold-forces"]
        if not threshold_stress:
            if "threshold-stress" in protocol_dict:
                threshold_stress = protocol_dict["threshold-stress"]
        if threshold_forces:
            parameters["md-max-force-tol"] = str(threshold_forces) + " eV/Ang"
        if threshold_stress:
            parameters["md-max-stress-tol"] = str(threshold_stress) + " eV/Ang**3"

        #Basis
        basis = self._get_basis(protocol, structure)

        #Kpoints
        if "kpoints" in protocol_dict:
            kpoints_mesh = KpointsData()
            kpoints_mesh.set_cell_from_structure(structure)
            kp_dict = protocol_dict["kpoints"]
            if "offset" in kp_dict:
                kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"], offset=kp_dict["offset"])
            else:
                kpoints_mesh.set_kpoints_mesh_from_density(distance=kp_dict["distance"])

        #Pseudo fam
        pseudo_fam = protocol_dict["pseudo_family"]

        #builder construction
        builder = SiestaBaseWorkChain.get_builder()
        builder.structure = structure
        builder.basis = Dict(dict=basis)
        builder.parameters = Dict(dict=parameters)
        if "kpoints" in protocol_dict:
            builder.kpoints = kpoints_mesh
        builder.pseudo_family = Str(pseudo_fam)
        builder.options = Dict(dict=calc_engines['relaxation']["options"])
        builder.code = load_code(calc_engines['relaxation']["code"])

        return builder
