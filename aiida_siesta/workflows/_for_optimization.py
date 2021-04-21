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
            if "$0" in str(value.get_dict()[key]):
                has_label = True
        if not has_label:
            return "the basis dictionary must contain the variable for optimizaton; they are '$'+int."


@calcfunction
def extract(out_par):
    """
    Extract the basis_enthalpy from the output of a SiestaBaseWorkChain if
    a Dict is passes. If a Str is passed instead, this means that the
    SiestaBaseWorkChain was not run and the workchain should return
    a big energy (we decided for 0.0).
    """
    if isinstance(out_par, orm.Str):
        return orm.Float(0.0)
    return orm.Float(out_par["basis_enthalpy"])


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
        #Do something to impose basis dict is present and the conyains the $
        spec.input('list_of_values', valid_type=orm.List)
        spec.input_namespace("simplex")
        spec.input('simplex.upper_boundaries', valid_type=orm.List)
        spec.input('simplex.lower_boundaries', valid_type=orm.List)
        spec.output("ene", valid_type=orm.Float)

        spec.inputs["siesta_base"]["basis"].validator = validator_basis_opt

        spec.outline(if_(cls.should_run_wc)(
            cls.prepare_inputs,
            cls.run_siesta_wc,
        ), cls.run_results)
        spec.exit_code(200, 'ERROR_WC', message='The SiestaBaseWorkChain failed')

    def should_run_wc(self):
        """
        Decides whether to run the SiestaBaseWorkChain. In fact, if the selected values for the
        variables are outside a range defined in input (see `simplex.upper/lower_boundaries`) we do
        not run siesta but return a high energy.
        """
        the_list = self.inputs.list_of_values.get_list()
        uppers = self.inputs.simplex.upper_boundaries.get_list()
        lowers = self.inputs.simplex.lower_boundaries.get_list()
        self.ctx.should_out = True

        for index, num in enumerate(the_list):
            if num > uppers[index] or num < lowers[index]:
                self.ctx.should_out = False
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
            for index, num in enumerate(self.inputs.list_of_values.get_list()):
                basis_dict[key] = basis_dict[key].replace("$" + str(index), str(num))

        self.ctx.inputs = inputs
        self.ctx.inputs['basis'] = orm.Dict(dict=basis_dict)

    def run_siesta_wc(self):
        """
        Run the SiestaBaseWorkChain.
        """

        self.report('Variables to optimiza correctly introduced in the basis dict.')

        running = self.submit(SiestaBaseWorkChain, **self.ctx.inputs)
        self.report(f'Launched SiestaBaseWorkChain<{running.pk}> to perform the siesta calculation.')

        return ToContext(workchain_base=running)

    def run_results(self):
        """
        Return the basis enthalpy in a single output node.
        """
        if self.ctx.should_out:
            if not self.ctx.workchain_base.is_finished_ok:
                extract_ene = extract(orm.Str("none"))
                self.out("ene", extract_ene)
                return self.exit_codes.ERROR_WC
            extract_ene = extract(self.ctx.workchain_base.outputs["output_parameters"])
        else:
            extract_ene = extract(orm.Str("none"))

        self.out("ene", extract_ene)

        return None
