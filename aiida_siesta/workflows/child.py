# -*- coding: utf-8 -*-

from aiida.work.workchain import while_, if_

from base import SiestaBaseWorkChain


class SiestaWorkChain(SiestaBaseWorkChain):
    @classmethod
    def create_outline(cls):
        outline = (
            cls.setup,
            cls.validate_pseudo_potentials,
            cls.should_setup,
            while_(cls.marry)(
                # NOTE rename base methods
                while_(cls.should_run_siesta)(
                    cls.run_siesta,
                    cls.inspect_siesta,
                ),
                # TODO decorate
                cls.decrement,
                cls.should_reset,
                #
                # cls.run_results,
            ),
            cls.run_results,
        )

        return(outline)

    def should_setup(self):
        self.ctx.sorry = 4

    def decrement(self):
        self.ctx.sorry -= 1

    def should_reset(self):
        self.ctx.is_finished = False
        self.ctx.iteration = 0

    def marry(self):
        if  self.ctx.sorry > 0:
            print 'Not Sorry'
            return True

        print 'Sorry'
        return False
