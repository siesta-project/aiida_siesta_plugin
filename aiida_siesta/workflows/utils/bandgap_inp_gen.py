from aiida_siesta.workflows.utils.base_inp_gen import BaseWorkChainInputsGenerator


class BandgapWorkChainInputsGenerator(BaseWorkChainInputsGenerator):

    #pylint: disable=signature-differs
    def get_builder(self, structure, calc_engines, protocol, path_generator, relaxation_type=None):

        return super().get_builder(structure, calc_engines, protocol, path_generator, relaxation_type)
