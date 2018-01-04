# -*- coding: utf-8 -*-

from aiida.work.workchain import while_, if_

from base import SiestaBaseWorkChain


class SiestaWorkChain(SiestaBaseWorkChain):
    @classmethod
    def create_outline(cls):
        outline = (
            cls._initial_setup,
            cls._validate_pseudo_potentials,
            cls.should_setup,
            while_(cls.ready)(
                *cls.run_scf(cls.decrement)
            ),
            cls._scf_results
        )

        return(outline)

    @classmethod
    def run_scf(cls, func):
        sequence = (cls._scf_reset,
                    while_(cls._should_run_scf)(cls._run_scf_cycle,
                                                cls._inspect_scf_cycle,),
                    func,)
        return sequence

    def should_setup(self):
        self.ctx.ready = 4

    def decrement(self):
        self.ctx.ready -= 1

    def ready(self):
        if  self.ctx.ready > 0:
            print 'Not Ready'
            return True

        print 'Ready'
        return False
