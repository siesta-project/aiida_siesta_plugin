# -*- coding: utf-8 -*-

from aiida.work.workchain import while_, if_

from base import SiestaBaseWorkChain


class SiestaBandsWorkChain(SiestaBaseWorkChain):
    @classmethod
    def create_outline(cls):
        outline = (
            cls._initial_setup,
            cls._validate_pseudo_potentials,
            cls._scf_reset,
            while_(cls._should_run_scf)(cls._run_scf_cycle,
                                        cls._inspect_scf_cycle,),
            cls.calculate_bands,
            cls.extra_results,
        )
        return(outline)

    def calculate_bands(self):
        if self.ctx.last_calc:
            local_inputs['parent_folder'] = self.ctx.last_calc.out.remote_folder

    def extra_results(self):
        self.out("bands_array", self.ctx.last_calc.out.bands_array)
