from aiida import orm
from aiida.engine import WorkChain, while_, ToContext
from aiida.common import AttributeDict
from aiida_siesta.workflows.simplex_basis import SimplexBasisOptimization


class TwoStepsBasisOpt(WorkChain):
    """
    Optimization that is more similar to the simplex code in the siesta utils.
    The optimization has two levels, a "marcrostep" that consists in the restart of
    a simplex with gradual reduction of the dimention of the initial simplex.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(
            SimplexBasisOptimization, exclude=('metadata', 'simplex.initial_step_fraction', 'simplex.full_simplex')
        )
        spec.input('macrostep.initial_lambda', valid_type=orm.Float, default=lambda: orm.Float(0.4))
        spec.input('macrostep.lambda_scaling_factor', valid_type=orm.Float, default=lambda: orm.Float(0.5))
        spec.input('macrostep.minimum_lambda', valid_type=orm.Float, default=lambda: orm.Float(0.01))
        spec.expose_outputs(SimplexBasisOptimization)
        spec.outline(
            cls.preprocess,
            while_(cls.should_run)(
                cls.run_simplex,
                cls.update_attributes,
            ),
            cls.run_results,
        )

    def preprocess(self):
        """
        Set some internal attributes.
        """
        self.ctx.current_lambda = self.inputs.macrostep.initial_lambda.value
        self.ctx.vars = self.inputs.simplex.initial_variables

    def should_run(self):
        """
        Run until the lambda is smaller than the final lambda.
        """
        return self.ctx.current_lambda > self.inputs.macrostep.minimum_lambda.value

    def run_simplex(self):
        """
        Run the simplex with the current lambda as `initial_step_fraction`.
        """
        min_lambda = self.inputs.macrostep.minimum_lambda.value
        self.report(f"Start macrostep with lambda {self.ctx.current_lambda}, minimum lambda {min_lambda}")
        siesta_base = self.inputs.siesta_base
        simpl = AttributeDict(self.inputs.simplex)
        simpl["initial_step_fraction"] = orm.Float(self.ctx.current_lambda)
        simpl["initial_variables"] = self.ctx.vars
        inputs = {"simplex": simpl, "siesta_base": siesta_base}
        run = self.submit(SimplexBasisOptimization, **inputs)
        return ToContext(simplex_wc=run)

    def update_attributes(self):
        """
        Updtae the value of lambda and the initial values of the variable for the next simplex.
        Please not we get the initial values of the variable from the `last_simplex` output
        of the last run optimization. This output is the only one that is returned when the
        optimization does not conclude in the max number of steps and it has the "best so far"
        set of values in the first position of the list.
        """
        self.ctx.current_lambda = self.ctx.current_lambda * self.inputs.macrostep.lambda_scaling_factor.value
        self.ctx.vars = orm.List(list=self.ctx.simplex_wc.outputs.last_simplex.get_list()[0])

    def run_results(self):
        """
        Return the outputs of the final simplex.
        """
        self.report("Lambda reached the minimum. End of two-steps optimization")
        self.out_many(self.exposed_outputs(self.ctx.simplex_wc, SimplexBasisOptimization))
