from aiida import orm
from aiida.engine import WorkChain, ToContext
from aiida_optimize.engines import NelderMead
from aiida_optimize import OptimizationWorkChain
from aiida_siesta.workflows._for_optimization import ForBasisOptWorkChain


class SimplexBasisOptimization(WorkChain):
    """
    Workchain running a simple NelderMead optimization (simplex) varing variables
    defined in the basis dictionaries.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(ForBasisOptWorkChain, exclude=('metadata', 'list_of_values'))
        spec.input('simplex.initial_variables', valid_type=orm.List)
        spec.input('simplex.max_iters', valid_type=orm.Int, default=lambda: orm.Int(40))
        spec.input('simplex.tolerance_function', valid_type=orm.Float, default=lambda: orm.Float(0.01))
        spec.input('simplex.initial_step_fraction', valid_type=orm.Float, default=lambda: orm.Float(0.4))
        spec.input('simplex.full_simplex', valid_type=orm.List, required=False)
        spec.output('last_simplex', valid_type=orm.List, required=True)
        spec.expose_outputs(OptimizationWorkChain)
        spec.outline(
            cls.preprocess,
            cls.run_optimizer,
            cls.outpts,
        )
        spec.exit_code(200, 'ERROR_OPT', message='The OptimizationWorkChain did not finished')
        #spec.exit_code(201, 'ERROR_FINAL_WC', message='The SiestaBaseWorkChain to obtain the bands failed')

    def preprocess(self):
        """
        In the preprocess, we create the initaial simplex.
        """
        import numpy as np

        simplex_inps = self.inputs.simplex
        if "full_simplex" in simplex_inps:
            self.ctx.simplex = simplex_inps["full_simplex"]
        else:
            simplex = []
            init_var_list = simplex_inps.initial_variables.get_list()
            ranges = np.array(simplex_inps.upper_boundaries) - np.array(simplex_inps.lower_boundaries)
            simplex.append(init_var_list)
            for index, num in enumerate(init_var_list):
                val = num + ranges[index] * simplex_inps.initial_step_fraction.value
                new_point = init_var_list.copy()
                new_point[index] = val
                simplex.append(new_point)
            self.ctx.simplex = simplex

        self.report("Successful set up of the initial simplex hypertetrahedron.")

    def run_optimizer(self):
        """
        Run the NelderMead optimization through the OptimizationWorkChain.
        """

        inputs = self.exposed_inputs(ForBasisOptWorkChain)

        siesta_base = inputs.siesta_base
        simpl = inputs.simplex

        result_opt = self.submit(
            OptimizationWorkChain,
            engine=NelderMead,
            engine_kwargs=orm.Dict(
                dict=dict(
                    simplex=self.ctx.simplex,
                    input_key="list_of_values",
                    result_key='ene',
                    xtol=10000,
                    ftol=self.inputs.simplex.tolerance_function.value,
                    max_iter=self.inputs.simplex.max_iters.value
                )
            ),
            evaluate_process=ForBasisOptWorkChain,
            evaluate={
                "simplex": {
                    "upper_boundaries": simpl.upper_boundaries,
                    "lower_boundaries": simpl.lower_boundaries
                },
                "siesta_base": {
                    **siesta_base
                }
            }
        )

        self.report("Launched ``OptimizationWorkChain`` with ``NelderMead`` engine")

        return ToContext(optimization=result_opt)

    def outpts(self):
        """
        Returning the outputs.
        """
        if not self.ctx.optimization.is_finished_ok:
            self.report("SimplexBasisOptimization is not finished succesfully, returning `last_simplex`.")
            self.out("last_simplex", self.ctx.optimization.outputs["engine_outputs__last_simplex"])
            return self.exit_codes.ERROR_OPT

        self.report("Concluded SimplexBasisOptimization succesfully, returning outputs.")

        self.out("last_simplex", self.ctx.optimization.outputs["engine_outputs__last_simplex"])

        self.out_many(self.exposed_outputs(self.ctx.optimization, OptimizationWorkChain))
