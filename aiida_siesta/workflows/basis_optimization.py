from aiida.engine import calcfunction, WorkChain, ToContext
from aiida import orm
from aiida_siesta.workflows.simplex_basis import SimplexBasisOptimization
from aiida_siesta.workflows.iterate import SiestaIterator


@calcfunction
def extract_pao_block(**ions):
    """
    From the outputs ions, create a basis block with the variables for the optimization.
    """
    ang_to_bohr = 1.8897161646321
    variables = {}
    card = "\n"
    for ion_name, ion in ions.items():  #pylint: disable=too-many-nested-blocks
        pao_mod = ion.get_pao_modifier()
        for n in pao_mod._gen_dict:  #pylint: disable=invalid-name
            for l in pao_mod._gen_dict[n]:  # noqa
                save_z1 = pao_mod._gen_dict[n][l][1]
                for z in pao_mod._gen_dict[n][l]:
                    string = ion_name + str(n) + str(l) + str(z)
                    if z == 1:
                        variables[string] = [1.5, 13.5, round(pao_mod._gen_dict[n][l][z] * ang_to_bohr, 6)]
                        pao_mod._gen_dict[n][l][z] = "$" + string
                    else:
                        init_fac = round(pao_mod._gen_dict[n][l][z] / save_z1, 6)
                        if init_fac >= 0.95:
                            init_fac = 0.94
                        variables[string] = [0.15, 0.95, init_fac]
                        pao_mod._gen_dict[n][l][z] = "-$" + string
        card += pao_mod.get_pao_block() + "\n"

    card += "%endblock pao-basis"
    basis = orm.Dict(
        dict={
            "%block pao-basis": card,
            "reparametrize-pseudos": True,
            "restricted-radial-grid": False,
            "pao-split-tail-norm": True
        }
    )

    simpl = orm.Dict(dict=variables)

    return {"basis": basis, "simplex": simpl}


class BasisOptimizationWorkChain(WorkChain):
    """
    WorkChain for basis optimization
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(
            SimplexBasisOptimization, exclude=('metadata', 'simplex.variables_dict', 'siesta_base.basis')
        )
        spec.input('basis_sizes', valid_type=orm.List, default=lambda: orm.List(list=["DZ", "DZP", "TZ"]))
        spec.input('sizes_monitored_quantity', valid_type=orm.Str, default=lambda: orm.Str("as_simplex"))
        spec.input('non_perturbative_pol', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.output('optimal_basis_block', valid_type=orm.Dict)
        spec.outline(
            cls.run_sizes,
            cls.run_optimizer,
            cls.run_results,
        )
        spec.exit_code(200, 'ERROR_SIZES_RUN', message='The SiestaIterator running pao sizes failed')
        spec.exit_code(201, 'ERROR_OPT', message='The basis optimization failed, probably not sufficient steps')

    def run_sizes(self):
        """
        Run the iterator trying different basis sizes
        """
        siesta_inputs = self.exposed_inputs(SimplexBasisOptimization).siesta_base
        if self.inputs.non_perturbative_pol.value:
            siesta_inputs["basis"] = orm.Dict(dict={"pao-non-perturbative-polarization-orbitals": True})
        sizes_list = self.inputs.basis_sizes.get_list()
        run_sizes = self.submit(
            SiestaIterator,
            **siesta_inputs,
            iterate_over={'pao-basis-size': sizes_list},
            batch_size=orm.Int(len(sizes_list))
        )
        self.report(f'Launched SiestaIterator<{run_sizes.pk}> to test basis sizes.')
        return ToContext(sizes_run=run_sizes)

    def run_optimizer(self):
        """
        Extract basis size with lowest "sizes_monitored_quantity", create a starting pao block and
        run the simplex optimizer varying all the cutoff radius.
        The "sizes_monitored_quantity" is the quantity used to discriminate which is the best
        basis size. A check is than on what basis size gives the smallest "sizes_monitored_quantity".
        This quantity can be defined in input. By default it is the one used for the simplex and
        specified in the input "simplex.output_name".
        """
        if not self.ctx.sizes_run.is_finished_ok:
            return self.exit_codes.ERROR_SIZES_RUN

        self.report('Concluded calculations to test basis sizes. Analyzing results.')

        if self.inputs.sizes_monitored_quantity.value == "as_simplex":
            monitoring_quantity = self.exposed_inputs(SimplexBasisOptimization).simplex.output_name.value
        else:
            monitoring_quantity = self.inputs.sizes_monitored_quantity.value
        monitoring_quantity_value = 10000000
        good_calc = None
        for calc in self.ctx.sizes_run.called:
            if calc.outputs.output_parameters[monitoring_quantity] <= monitoring_quantity_value:
                monitoring_quantity_value = calc.outputs.output_parameters[monitoring_quantity]
                good_calc = calc

        self.report(f'Looking at {monitoring_quantity}, {good_calc.inputs.basis["paobasissize"]} is favourable')

        if "ion_files" in good_calc.get_outgoing().nested():
            ions = good_calc.get_outgoing().nested()["ion_files"]

        pao_opt_all = extract_pao_block(**ions)

        siesta_inputs = self.exposed_inputs(SimplexBasisOptimization).siesta_base
        siesta_inputs["basis"] = pao_opt_all["basis"]
        simplex_inputs = self.exposed_inputs(SimplexBasisOptimization).simplex
        simplex_inputs["variables_dict"] = pao_opt_all["simplex"]
        simp_wc_inputs = {
            'siesta_base': siesta_inputs,
            'simplex': simplex_inputs,
        }
        run_opt = self.submit(SimplexBasisOptimization, **simp_wc_inputs)
        self.report(f'Launched SimplexBasisOptimization<{run_opt.pk}>.')
        return ToContext(simplex_run=run_opt)

    def run_results(self):
        """
        Extract outputs
        """
        if not self.ctx.simplex_run.is_finished_ok:
            return self.exit_codes.ERROR_OPT

        self.report("Concluded SimplexBasisOptimization succesfully, returning outputs.")

        opt_for_basis_wc = orm.load_node(self.ctx.simplex_run.outputs.optimal_process_uuid.value)

        opt_basis = opt_for_basis_wc.get_outgoing(node_class=orm.WorkChainNode).one().node.inputs.basis

        self.out("optimal_basis_block", opt_basis)
