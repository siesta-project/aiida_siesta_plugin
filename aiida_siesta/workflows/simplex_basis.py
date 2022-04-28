from aiida import orm
from aiida.engine import WorkChain, ToContext
from aiida_optimize.engines import NelderMead
from aiida_optimize import OptimizationWorkChain
from aiida_siesta.workflows._for_optimization import ForBasisOptWorkChain


def validate_variables_dict(value, _):
    """
    Check that each key in vriables_dict respect the format.
    """
    if value:
        dime = len(value.get_dict().keys())
        for k, v in value.get_dict().items():
            if not isinstance(v, (tuple, list)):
                return "the values for each key must be list or tuple."
            if len(v) not in [2, 3, dime + 1 + 2]:
                messag = (
                    f"for `{k}` of the `variables_dict` input. The list/tuple must have lenght 2, 3 or {dime+1+2}. " +
                    f"First number is the minimum for `{k}`, the second number the maximum, the third (optional) " +
                    "the initial value, the others (optional) the numbers to create the initial simplex."
                )
                return messag
            if not all([isinstance(el, (float, int)) for el in v]):
                return f"the values for each key {k} must be list/tuple of floats or integers."


class SimplexBasisOptimization(WorkChain):
    """
    Workchain running a simple NelderMead optimization (simplex) varing variables
    defined in the basis dictionaries.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(
            ForBasisOptWorkChain,
            exclude=('metadata', 'the_values', 'the_names', 'upper_bounds', 'lower_bounds', 'out_name')
        )
        spec.input_namespace("simplex")
        spec.input('simplex.variables_dict', valid_type=orm.Dict, validator=validate_variables_dict)
        spec.input('simplex.max_iters', valid_type=orm.Int, default=lambda: orm.Int(40))
        spec.input('simplex.tolerance_function', valid_type=orm.Float, default=lambda: orm.Float(0.01))
        spec.input('simplex.initial_step_fraction', valid_type=orm.Float, default=lambda: orm.Float(0.4))
        spec.input('simplex.output_name', valid_type=orm.Str, default=lambda: orm.Str("basis_enthalpy"))
        #spec.input('simplex.full_simplex', valid_type=orm.List, required=False)
        spec.output('last_simplex', valid_type=orm.List, required=True)
        spec.expose_outputs(OptimizationWorkChain, exclude="engine_outputs")
        spec.outline(
            cls.preprocess,
            cls.run_optimizer,
            cls.outpts,
        )
        spec.exit_code(200, 'ERROR_OPT', message='The OptimizationWorkChain did not finished')
        #spec.exit_code(201, 'ERROR_FINAL_WC', message='The SiestaBaseWorkChain to obtain the bands failed')

    def preprocess(self):
        """
        In the preprocess, we transform the `variables_dict` info into lists. This
        is necessary to exploit the OptimizationWorkChain features. It accepts the variables
        values as list in a separate input. Also the initaial simplex is created since it is
        explicitly needed by OptimizationWorkChain.
        """
        import numpy as np
        import random

        simplex_inps = self.inputs.simplex
        dime = len(simplex_inps.variables_dict.get_dict().keys())
        self.ctx.names = []
        self.ctx.lower = []
        self.ctx.upper = []
        init = []
        others = []
        for ind in range(dime):
            others.append([])
        for k, v in simplex_inps.variables_dict.get_dict().items():
            self.ctx.names.append(k)
            self.ctx.lower.append(v[0])
            self.ctx.upper.append(v[1])
            if len(v) > 2:
                init.append(v[2])
            else:
                init.append(round(random.uniform(v[0], v[1]), 5))
            if len(v) > 3:
                for ind in range(dime):
                    others[ind].append(v[ind + 2])
        if all([len(others[ind]) == dime for ind in range(dime)]):
            #Means comlete simplex is passed
            self.ctx.simplex = others.insert(0, init)
        else:
            simplex = []
            ranges = np.array(self.ctx.upper) - np.array(self.ctx.lower)
            simplex.append(init)
            for index, num in enumerate(init):
                if num > self.ctx.upper[index] / 2 + 0.0001:
                    val = num - ranges[index] * simplex_inps.initial_step_fraction.value
                else:
                    val = num + ranges[index] * simplex_inps.initial_step_fraction.value
                new_point = init.copy()
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

        result_opt = self.submit(
            OptimizationWorkChain,
            engine=NelderMead,
            engine_kwargs=orm.Dict(
                dict=dict(
                    simplex=self.ctx.simplex,
                    input_key="the_values",
                    result_key='ene',
                    xtol=10000,
                    ftol=self.inputs.simplex.tolerance_function.value,
                    max_iter=self.inputs.simplex.max_iters.value
                )
            ),
            evaluate_process=ForBasisOptWorkChain,
            evaluate={
                "upper_bounds": orm.List(list=self.ctx.upper),
                "lower_bounds": orm.List(list=self.ctx.lower),
                "the_names": orm.List(list=self.ctx.names),
                "out_name": self.inputs.simplex.output_name,
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
