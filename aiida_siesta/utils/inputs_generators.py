from .generator_absclass import InputsGenerator


class SiestaCalculationInputsGenerator(InputsGenerator):
    """
    This class has two main purposes:
    1) Provide a method (get_filled_builder) that returns a builder for SiestaCalculation
       with pre-compiled inputs according to a protocol and some relaxation/bands/spin options.
    2) Implement few methods that can be used to explore the relaxation types available,
       the bands options and the spin options (protocols options are inherited from ProtocolManager)
    When instanciated, it requires the SiestaCalculation class as input. This is needed to avoid
    cyclic dependence.
    """

    _calc_types = {
        "siesta": {
            'code': 'Put here the code name, must be for plugin siesta.siesta',
            'options': 'Put here the computational options for running the relaxation, following the usual '
            'aiida schema'
        }
    }

    _relax_types = {
        "atoms_only": "The lattice shape and volume are fixed, only the atomic positions are relaxed",
        "variable_cell": "The lattice is relaxed together with the atomic coordinates",
        "constant_volume": "The cell volume is kept constant in a variable-cell relaxation"
    }
    _bands_path_generators = {
        "seekpath":
        "Cell converted always in the primitive cell. More relaiable",
        "legacy":
        "To be used only if it is important to keep unchanged the cell shape"
        "of the input sturcture (no primitive cell); not guaranteed to give correct paths",
    }
    _spins = {
        "polarized": "Collinear spin option",
        "non-collinear": "Non-collinear spin option",
        "spin-orbit": "Spin-orbit is activated",
    }

    #Some methods to return info about options that the get_filled_builder of this class
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

    def get_spins(self):
        return list(self._spins.keys())

    def get_spin_info(self, key):
        try:
            return self._spins[key]
        except KeyError:
            raise ValueError("Wrong spin type: no spin with name {} implemented".format(key))

    #pylint: disable=too-many-statements,arguments-differ
    def get_inputs_dict(  # noqa: MC0001  - is mccabe too complex funct -
        self, structure, calc_engines, protocol, bands_path_generator=None, relaxation_type=None, spin=None
    ):
        """
        Method return a dictionary with the inputs of a SiestaBaseWorkChain, obtained according
        to the protocol, structure, calc_engines passed by the user (optionally bands_path_generator
        relaxation_type and spin as well)
        """

        from aiida.orm import (Dict, KpointsData)
        from aiida.orm import load_code
        from aiida.tools import get_explicit_kpoints_path

        #Checks
        if spin is not None:
            if spin not in self.get_spins():
                raise ValueError("No `spin` with name {} implemented".format(bands_path_generator))
        if bands_path_generator is not None:
            if bands_path_generator not in self.get_bands_path_generators():
                raise ValueError("No `bands_path_generator` with name {} implemented".format(bands_path_generator))
        if relaxation_type is not None:
            if relaxation_type not in self.get_relaxation_types():
                raise ValueError("No `relaxation_type` with name {} implemented".format(relaxation_type))
        if not self.is_valid_protocol(protocol):
            import warnings
            defpro = self.get_default_protocol_name()
            warnings.warn("no protocol implemented with name `{0}`, using default `{1}`".format(protocol, defpro))
            protocol = defpro
        if 'siesta' not in calc_engines:
            raise ValueError("Wrong syntax in `calc_engines`. Check method `how_to_pass_computation_options`.")

        #Bandskpoints (might change structure)
        if bands_path_generator is not None:
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
        if spin is not None:
            sp_parameters = self._add_spin_options(protocol, prot_param)
            sp_parameters["spin"] = spin
        else:
            sp_parameters = prot_param.copy()
        #Add relaxation options in case requested
        if relaxation_type is not None:
            parameters = self._add_relaxation_options(protocol, sp_parameters)
            parameters["md-type-of-run"] = "cg"
            parameters["md-num-cg-steps"] = 100
            if relaxation_type == "variable_cell":
                parameters["md-variable-cell"] = True
            if relaxation_type == "constant_volume":
                parameters["md-variable-cell"] = True
                parameters["md-constant-volume"] = True
        else:
            parameters = sp_parameters.copy()

        #Basis
        basis = self._get_basis(protocol, ok_structure)

        #Kpoints (might not be present, for molecules for instance)
        kpoints_mesh = self._get_kpoints(protocol, ok_structure)

        #Pseudo fam
        pseudos = self._get_pseudos(protocol, ok_structure)

        #Computational resources
        options = calc_engines['siesta']["options"]
        code = load_code(calc_engines['siesta']["code"])

        inputs = {
            'structure': ok_structure,
            'parameters': Dict(dict=parameters),
            'code': code,
            'basis': Dict(dict=basis),
            'pseudos': pseudos,
            'metadata': {
                "options": options
            },
        }

        #bandskpoints and kpoints are optional, possible return None
        if kpoints_mesh is not None:
            inputs['kpoints'] = kpoints_mesh
        if bandskpoints is not None:
            inputs['bandskpoints'] = bandskpoints

        return inputs

    # pylint: disable=arguments-differ
    def get_filled_builder(
        self, structure, calc_engines, protocol, bands_path_generator=None, relaxation_type=None, spin=None
    ):
        """
        Return a builder for the SiestaBaseWorkChain, pre filled according to the protocol,
        structure, calc_engines passed by the user (optionally bands_path_generator
        relaxation_type and spin as well). Ready to be submitted.
        """

        inp_dict = self.get_inputs_dict(structure, calc_engines, protocol, bands_path_generator, relaxation_type, spin)

        builder = self._fill_builder(inp_dict)

        return builder


class BaseWorkChainInputsGenerator(SiestaCalculationInputsGenerator):
    """
    Inputs generator for the SiestaBaseWorkChain, it must receive
    a SiestaBaseWorkChain as argument when instanciated.
    Makes use of the methods of the SiestaCalculationInputsGenerator, only the
    handling of the options needs to be changed.
    """

    def get_inputs_dict(
        self, structure, calc_engines, protocol, bands_path_generator=None, relaxation_type=None, spin=None
    ):

        from aiida.orm import Dict

        inps = super().get_inputs_dict(structure, calc_engines, protocol, bands_path_generator, relaxation_type, spin)

        inps['options'] = Dict(dict=inps['metadata']['options'])
        inps['metadata'] = {}

        return inps


#class BandgapWorkChainInputsGenerator(BaseWorkChainInputsGenerator):
#    """
#    Inputs generator for the BandgapWorkChain, makes use of the methods
#    of the BaseWorkChainInputsGenerator, only the __init__ requires the correct
#    workchain class in input and the `get_inputs_dict` implements a check for the
#    presence of `bands_path_generator`. In fact the band calculation is required.
#    """
#
#    def get_inputs_dict(
#        self, structure, calc_engines, protocol, bands_path_generator=None, relaxation_type=None, spin=None
#    ):
#
#        if not bands_path_generator:
#            raise RuntimeError(
#                'Method `get_inputs_dict` of class `{0}` requires `bands_path_generator`'.format(
#                    self.__class__.__name__
#                )
#            )
#
#        return super().get_inputs_dict(structure, calc_engines, protocol, bands_path_generator, relaxation_type, spin)


class EosWorkChainInputsGenerator(BaseWorkChainInputsGenerator):
    """
    Inputs generator for the FixedCellEoS WorkChain, makes use of the methods
    of the BaseWorkChainInputsGenerator, only the __init__ requires the correct
    workchain class in input and the `get_inputs_dict` implements a check for the
    presence of an unsupported relaxation.
    """

    _relax_types = {
        "atoms_only": "The lattice shape and volume are fixed, only the atomic positions are relaxed",
        "constant_volume": "The cell volume is kept constant in a variable-cell relaxation"
    }


class StmWorkChainInputsGenerator(BaseWorkChainInputsGenerator):
    """
    Inputs generator for the STMWorkChain, makes use of the methods
    of the BaseWorkChainInputsGenerator, but, in addition to the __init__ requiring the correct
    workchain class in input, also the _calc_types needs to be modified to have the possibility
    to pass computational resources for the stm plugin. The `get_inputs_dict` implements the selection
    of many more inputs related to the STM. Same the `get_filled_builder`.
    """

    _calc_types = {
        "siesta": {
            'code': 'Put here the code name, must be for plugin siesta.siesta',
            'options': 'Put here the computational options for running the siesta code, following the usual '
            'aiida schema'
        },
        "stm": {
            'code': 'Put here the code name, must be for plugin siesta.stm',
            'options': 'Put here the computational options for running the stm, following the usual aiida schema'
        }
    }
    _stm_modes = {
        "constant-height":
        "primitive cell always used. Always ok",
        "constant-current":
        "To be used only if it is important to keep unchanged the cell shape"
        "of the input sturcture (no primitive cell). Not guaranteed to give correct paths",
    }

    def get_stm_modes(self):
        return list(self._stm_modes.keys())

    def get_stm_mode_info(self, key):
        try:
            return self._stm_modes[key]
        except KeyError:
            raise ValueError("Wrong stm mode: no stm_mode with name {} implemented".format(key))

    def get_stm_value_info(self):  # pylint: disable=no-self-use
        return "Value of height in Ang or value of current in e/bohr**3 (float)"

    # pylint: disable=arguments-differ
    def get_inputs_dict(
        self,
        structure,
        calc_engines,
        protocol,
        stm_mode,
        stm_value,
        bands_path_generator=None,
        relaxation_type=None,
        spin=None
    ):

        from aiida.orm import (Dict, Float, Str)
        from aiida.orm import load_code

        siesta_in = super().get_inputs_dict(
            structure, calc_engines, protocol, bands_path_generator, relaxation_type, spin
        )

        if stm_mode not in self._stm_modes:
            raise ValueError("Wrong stm mode: no stm_mode with name {} implemented".format(stm_mode))
        try:
            float(stm_value)
        except ValueError:
            raise RuntimeError("Wrong stm value: it must be a `float`")
        if 'stm' not in calc_engines:
            raise ValueError("Wrong syntax in `calc_engines`. Check method `how_to_pass_computation_options`.")

        stm_spin = "none"
        if spin:
            if spin == "polarized":
                stm_spin = "collinear"
            else:
                stm_spin = "non-collinear"

        siesta_in["stm_spin"] = Str(stm_spin)
        siesta_in["stm_mode"] = Str(stm_mode)
        siesta_in["stm_value"] = Float(stm_value)
        siesta_in["emin"] = Float(-6.5)
        siesta_in["emax"] = Float(+0.1)
        siesta_in["stm_options"] = Dict(dict=calc_engines['stm']["options"])
        siesta_in["stm_code"] = load_code(calc_engines['stm']["code"])

        return siesta_in

    # pylint: disable=arguments-differ
    def get_filled_builder(
        self,
        structure,
        calc_engines,
        protocol,
        stm_mode,
        stm_value,
        bands_path_generator=None,
        relaxation_type=None,
        spin=None
    ):

        inp_dict = self.get_inputs_dict(
            structure, calc_engines, protocol, stm_mode, stm_value, bands_path_generator, relaxation_type, spin
        )

        builder = self._fill_builder(inp_dict)

        return builder
