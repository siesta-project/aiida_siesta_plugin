from aiida_siesta.workflows.utils.base_inp_gen import BaseWorkChainInputsGenerator


class BandgapWorkChainInputsGenerator(BaseWorkChainInputsGenerator):

    def get_inputs_dict(self, structure, calc_engines, protocol, bands_path_generator=None, relaxation_type=None):

        if not bands_path_generator:
            RuntimeError("BandgapWorkChainInputsGenerator requires bands_path_generator")

        return super().get_inputs_dict(structure, calc_engines, protocol, bands_path_generator, relaxation_type)

    def get_builder(self, structure, calc_engines, protocol, bands_path_generator=None, relaxation_type=None):

        if not bands_path_generator:
            RuntimeError("BandgapWorkChainInputsGenerator requires bands_path_generator")

        return super().get_builder(structure, calc_engines, protocol, bands_path_generator, relaxation_type)
