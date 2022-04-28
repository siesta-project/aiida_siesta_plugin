from aiida import orm
from aiida.engine import WorkChain, calcfunction, ToContext, if_
from aiida.common import AttributeDict
from aiida_siesta.workflows.base import SiestaBaseWorkChain


def validator_basis_opt(value, _):
    """
    Integrate the `basis` validator of the SiestaBaseWorkChain with a check
    on the fact that we shuould have at least one variable in the basis dict,
    meaning a quantity that is varied in the optimization procedure.
    We decided to recognise the variables using the character '$' followed by
    an integer, starting from 0.
    """
    if value:
        orig_out = SiestaBaseWorkChain.spec().inputs["basis"].validator(value, _)
        if orig_out is not None:
            return orig_out
        has_label = False
        for key in value.get_dict():
            if "$" in str(value.get_dict()[key]):
                has_label = True
        if not has_label:
            return "the basis dictionary must contain the variable for optimizaton; they are '$'+str."


@calcfunction
def extract(out_par, out_name):
    """
    Extract the basis_enthalpy from the output of a SiestaBaseWorkChain if
    a Dict is passes. If a Str is passed instead, this means that the
    SiestaBaseWorkChain was not run and the workchain should return
    a big energy (we decided for 0.0).
    """
    if isinstance(out_par, orm.Str):
        return orm.Float(0.0)
    return orm.Float(out_par[out_name.value])


class ForBasisOptWorkChain(WorkChain):
    """
    Class wrapping the SiestaBaseWorkChain with the scope of
    assigning the variables of the optimization process to the right
    input of the SiestaBaseWorkChain.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(SiestaBaseWorkChain, namespace="siesta_base", exclude=('metadata',))
        spec.input('the_values', valid_type=orm.List)
        spec.input('the_names', valid_type=orm.List)
        spec.input('upper_bounds', valid_type=orm.List)
        spec.input('lower_bounds', valid_type=orm.List)
        spec.input('out_name', valid_type=orm.Str)
        spec.output("ene", valid_type=orm.Float)

        spec.inputs["siesta_base"]["basis"].validator = validator_basis_opt

        spec.outline(if_(cls.should_run_wc)(
            cls.prepare_inputs,
            cls.run_siesta_wc,
        ), cls.run_results)

    def should_run_wc(self):
        """
        Decides whether to run the SiestaBaseWorkChain. In fact, if the selected values for the
        variables are outside a range defined in input (see `simplex.upper/lower_boundaries`) we do
        not run siesta but return a high energy.
        """
        self.ctx.vals = self.inputs.the_values.get_list()
        self.ctx.names = self.inputs.the_names.get_list()
        uppers = self.inputs.upper_bounds.get_list()
        lowers = self.inputs.lower_bounds.get_list()

        for index, num in enumerate(self.ctx.vals):
            if num > uppers[index] or num < lowers[index]:
                return False

        return True

    def prepare_inputs(self):
        """
        Method that inserts the correct variable values inside the basis dict.
        Relies on the decition to have '$' character.
        """
        inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain, namespace='siesta_base'))

        basis_dict = inputs["basis"].get_dict()
        for key in basis_dict.keys():
            for index, num in enumerate(self.ctx.vals):
                if isinstance(basis_dict[key], str):  #Otherwise fail if number or bool
                    basis_dict[key] = basis_dict[key].replace("$" + str(self.ctx.names[index]), str(num))

        self.ctx.inputs = inputs
        self.ctx.inputs['basis'] = orm.Dict(dict=basis_dict)

    def run_siesta_wc(self):
        """
        Run the SiestaBaseWorkChain.
        """

        self.report('Variables to optimise correctly introduced in the basis dict.')

        running = self.submit(SiestaBaseWorkChain, **self.ctx.inputs)
        self.report(f'Launched SiestaBaseWorkChain<{running.pk}> to perform the siesta calculation.')

        return ToContext(workchain_base=running)

    def run_results(self):
        """
        Return the basis enthalpy in a single output node.
        """
        if self.should_run_wc():
            if not self.ctx.workchain_base.is_finished_ok:
                extract_ene = extract(orm.Str("none"), self.inputs.out_name)
            else:
                extract_ene = extract(self.ctx.workchain_base.outputs["output_parameters"], self.inputs.out_name)
        else:
            extract_ene = extract(orm.Str("none"), self.inputs.out_name)

        self.out("ene", extract_ene)
