# -*- coding: utf-8 -*-

from aiida.work.workchain import while_, if_

from base import SiestaBaseWorkChain


class SiestaWorkChain(SiestaBaseWorkChain):
    @classmethod
    def create_outline(cls):
        outline = (
            cls._initial_setup,
            cls._validate_pseudo_potentials,
            cls._scf_reset,
            while_(cls._should_run_scf)(cls._run_scf_cycle,
                                        cls._inspect_scf_cycle,),
            cls._scf_results,
        )
        return(outline)
