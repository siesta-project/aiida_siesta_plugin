from aiida.engine import WorkChain, ToContext, if_
from aiida import orm
from aiida_siesta.workflows.simplex_basis import SimplexBasisOptimization
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.workflows.iterate import SiestaIterator
from aiida_siesta.utils.tkdict import FDFDict


def set_variables_in_pao_mod(ion_name, variables, n, l, pao_mod, global_splitnorm, charge_conf):  # noqa
    """
    Changes the pao_mod dictionary to include variables for obtimization
    """
    ang_to_bohr = 1.8897161646321
    save_z1 = pao_mod._gen_dict[n][l][1]
    for z in pao_mod._gen_dict[n][l]:
        string = ion_name + str(n) + str(l) + str(z)
        if z == 1:
            suggest_radius = round(pao_mod._gen_dict[n][l][z] * ang_to_bohr, 6)
            if suggest_radius * 2 < 12.5:
                variables[string] = [1.5, 12.5, suggest_radius]
            else:
                variables[string] = [1.5, suggest_radius * 2, suggest_radius]
            pao_mod._gen_dict[n][l][z] = "$" + string
            if charge_conf:
                if pao_mod._gen_occu[n][l][z] == 0.0:
                    string2 = ion_name + str(n) + str(l) + "Zconf"
                    variables[string2] = [1, 6, 3]
                    pao_mod._conf_dict = {"Q": {n: {l: ["$" + string2]}}}
        else:
            if global_splitnorm:
                pao_mod._gen_dict[n][l][z] = 0.0
            else:
                init_fac = round(pao_mod._gen_dict[n][l][z] / save_z1, 6)
                if init_fac >= 0.95:
                    init_fac = 0.94
                variables[string] = [0.15, 0.95, init_fac]
                pao_mod._gen_dict[n][l][z] = "-$" + string

    if global_splitnorm:
        variables["splitnorm"] = [0.05, 0.65, 0.15]

    return variables


def extract_pao_block(global_splitnorm, global_enrgyshift, charge_conf, **ions):
    """
    From the outputs ions, create a basis block with the variables for the optimization.
    It also selects range of the variables and initial point. This choice is not so
    trivial. The most natural decision would be to start from the radii suggested by
    the basis size choice but not sure for simplex it is good to have an optimal
    value as a starting point. The range as well is not trivial since the nature of
    first orbital is different compared to others. After several attempts I decide to
    go with:
    * first zeta: start is suggested radius by sizes calculation (rs), min is 1.5, max is
      2*rs (except when 2*rs < 12.5, then max is 12.5)
    * second or more zetas: start is 1/z (so 0.5 for z=2, 0.333.. for z=3, ...), min is 0.15,
      max is 0.95
    """
    variables = {}
    card = "\n"
    for ion_name, ion in ions.items():
        pao_mod = ion.get_pao_modifier()
        for n in pao_mod._gen_dict:  #pylint: disable=invalid-name
            for l in pao_mod._gen_dict[n]:  # noqa
                if global_enrgyshift:
                    for z in pao_mod._gen_dict[n][l]:
                        pao_mod._gen_dict[n][l][z] = 0.0
                else:
                    variables = set_variables_in_pao_mod(
                        ion_name, variables, n, l, pao_mod, global_splitnorm, charge_conf
                    )
        card += pao_mod.get_pao_block() + "\n"

    card += "%endblock pao-basis"

    basis_dict = {"%block pao-basis": card}
    if global_enrgyshift:
        basis_dict.update({"pao-split-norm": "$splitnorm", "pao-energy-shift": "$energyshift meV"})
        variables["splitnorm"] = [0.05, 0.65, 0.15]
        variables["energyshift"] = [0.5, 50.5, 5]
    if global_splitnorm:
        basis_dict.update({"pao-split-norm": "$splitnorm"})

    #simpl = orm.Dict(dict=variables)
    return basis_dict, variables
    #return {"basis": basis, "simplex": simpl}


def add_orbitals(orb_dict, **ions):
    """
    From the outputs ions, create a basis block with the variables for the optimization.
    It also selects range of the variables and initial point.
    First zeta radia are variables and the split-norm value. Also a charge confinament is added
    for empty orbitals. The charge value is a variable to optimize.
    """
    l_dict = {"s": 0, "p": 1, "d": 2, "f": 3}

    def check_max_n_with_particular_l(dictionary, the_l):
        """
        Find the maximum n with a particular l
        """
        max_n = -1
        for the_n in dictionary:
            if the_l in dictionary[the_n]:
                if the_n > max_n:
                    max_n = the_n
        return max_n

    card = "\n"
    one_changed = False
    for ion_name, ion in ions.items():
        to_change = orb_dict.get(ion_name, None)
        if isinstance(to_change, str):
            to_change = [to_change]
        pao_mod = ion.get_pao_modifier()
        for orb in to_change:
            n = int(orb[0])  #pylint: disable=invalid-name
            l = orb[1]  # noqa
            z = int(orb[2])
            #In order to avoid the error "orbital not bound, need to set a radius explicitely"
            #We set the radius of the highest occupied l-1 shell. If not present l-2.
            ok_l = l_dict[l] - 1
            max_n = check_max_n_with_particular_l(pao_mod._gen_dict, ok_l)
            if max_n == -1:
                ok_l = l_dict[l] - 2
                max_n = check_max_n_with_particular_l(pao_mod._gen_dict, ok_l)
            if max_n == -1:
                rad = 0.0
            else:
                rad = pao_mod._gen_dict[max_n][ok_l][1]
            try:
                pao_mod.add_orbital("Ang", rad, n, l_dict[l])
                one_changed = True
            except ValueError:
                continue
            if z == 2:
                pao_mod.add_orbital("Ang", 0.0, n, l_dict[l], 2)

        card += pao_mod.get_pao_block() + "\n"

    card += "%endblock pao-basis"

    if not one_changed:
        return None

    return card


def validate_basis_opt(value, _):
    """
    Integrate the `basis` validator of the SiestaBaseWorkChain with a check
    on the presence of "%block pao-basis" and 'pao-basis-size'. They are not allowed because
    they are set automatically by the workflow.
    Also remove the check on the presence of a '$' character, that was implemented
    in the ForBasisOptWorkChain class.
    """
    if value:
        orig_out = SiestaBaseWorkChain.spec().inputs["basis"].validator(value, _)
        if orig_out is not None:
            return orig_out
        for key in value.get_dict():
            if FDFDict.translate_key(key) == "%block paobasis":
                return "the `basis` dictionary can not contain the key '%block pao-basis', is set by the algorithm"
            if FDFDict.translate_key(key) == "paobasissize":
                return "the `basis` dictionary can not contain the key 'pao-basis-size', is set by the algorithm"


def validate_opt_schema(value, _):
    """
    Validate the `optimization_schema` port
    """
    if value:
        if value["global_energy_shift"].value and not value["global_split_norm"].value:
            return "Not supported global_energy_shift=True AND global_split_norm=False. Change `optimization_schema`."
        if value["global_energy_shift"].value and value["charge_confinement"].value:
            return "charge_confinement incompatible with global_energy_shift=True. Change `optimization_schema`."


def validate_add_orbital(value, _):  #pylint: disable=too-many-return-statements
    """
    Validate the `add_orbital` input port.
    Imposes a signature of the dict of this kind:
        {"Ca": ["3d1","4f1"], "Sr": "4d1"}
    Meaning the key is an element and value is the orbitals to add.
    At the moment we have suppot only for elemental symbols? Should we extend to
    kind names?
    """
    from aiida.common.constants import elements

    message = "The value of each item in the `add_orbital` dict must be a list of strings 'nlz', e.g. ['3d1']."

    if value:
        symbols = [v['symbol'] for v in elements.values()]
        for key, val in value.get_dict().items():
            if key not in symbols:
                return "Every key of `add_orbital` dictionary must be a symbol of the periodic table"
            if isinstance(val, str):
                val = [val]
            if not isinstance(val, list):
                return message
            for orb in val:
                if not len(orb) == 3:
                    return message
                try:
                    int(orb[0])
                except ValueError:
                    return message
                try:
                    int(orb[2])
                except ValueError:
                    return message
                if orb[1] not in ["s", "p", "d", "f"]:
                    return "The 'l' in the string 'nlz' in `add_orbital` dict has s, p, d and f as allowed values."
                if int(orb[2]) not in [1, 2]:
                    return "The 'z' in the string 'nlz' in `add_orbital` dict  has 1 and 2 as allowed values."


class BasisOptimizationWorkChain(WorkChain):
    """
    WorkChain for basis optimization
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(SimplexBasisOptimization, exclude=('metadata', 'simplex.variables_dict'))
        spec.input('basis_sizes', valid_type=orm.List, default=lambda: orm.List(list=["DZ", "DZP", "TZ"]))
        spec.input('sizes_monitored_quantity', valid_type=orm.Str, default=lambda: orm.Str("as_simplex"))
        spec.input('add_orbital', valid_type=orm.Dict, required=False, validator=validate_add_orbital)
        spec.input_namespace('optimization_schema', validator=validate_opt_schema)
        spec.input('optimization_schema.global_energy_shift', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.input('optimization_schema.global_split_norm', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.input('optimization_schema.charge_confinement', valid_type=orm.Bool, default=lambda: orm.Bool(False))

        spec.inputs["siesta_base"]["basis"].validator = validate_basis_opt

        spec.output('optimal_basis_block', valid_type=orm.Dict)

        spec.outline(
            cls.run_sizes,
            cls.check_sizes,
            if_(cls.should_add_orbitals)(cls.run_extra, cls.check_extra),
            cls.run_optimizer,
            cls.run_results,
        )

        spec.exit_code(200, 'ERROR_SIZES_RUN', message='The SiestaIterator running pao sizes failed')
        spec.exit_code(201, 'ERROR_OPT', message='The basis optimization failed, probably not sufficient steps')
        spec.exit_code(202, 'ERROR_EXTRA_RUN', message='The basis optimization failed, probably not sufficient steps')

    def run_sizes(self):
        """
        Run the iterator trying different basis sizes
        """
        siesta_inputs = self.exposed_inputs(SimplexBasisOptimization).siesta_base
        sizes_list = self.inputs.basis_sizes.get_list()
        run_sizes = self.submit(
            SiestaIterator,
            **siesta_inputs,
            iterate_over={'pao-basis-size': sizes_list},
            batch_size=orm.Int(len(sizes_list))
        )
        self.report(f'Launched SiestaIterator<{run_sizes.pk}> to test basis sizes.')
        return ToContext(sizes_run=run_sizes)

    def check_sizes(self):
        """
        Check what is the size with minimum of `sizes_monitored_quantity`.
        """
        if not self.ctx.sizes_run.is_finished_ok:
            return self.exit_codes.ERROR_SIZES_RUN
        for calc in self.ctx.sizes_run.called:
            if not calc.is_finished_ok:
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

        self.ctx.monitoring_quantity = monitoring_quantity
        return ToContext(good_sizes_calc=good_calc)

    def should_add_orbitals(self):
        """
        Check on the need of an extra step to check extra orbitals
        """
        if 'add_orbital' not in self.inputs:
            return False

        struct = self.exposed_inputs(SimplexBasisOptimization).siesta_base.structure
        present_kinds = [kind.symbol for kind in struct.kinds]

        for key in self.inputs['add_orbital'].get_dict():
            if key in present_kinds:
                return True

        message = (
            'Port `add_orbital` is set in input, but the keys of its dictionary does ' +
            'not correnspond to any kind in input structure. Skip check on additional orbitals.'
        )

        self.report(message)

        return False

    def run_extra(self):
        """
        Run an extra calculation addind d orbitals
        """
        if "ion_files" in self.ctx.good_sizes_calc.get_outgoing().nested():
            ions = self.ctx.good_sizes_calc.get_outgoing().nested()["ion_files"]

        pao_block = add_orbitals(self.inputs.add_orbital.get_dict(), **ions)

        if pao_block is None:
            self.report("Input `add_orbital` included, but no orbital could be succesfully added. Skip")
            return None

        siesta_inputs = self.exposed_inputs(SimplexBasisOptimization).siesta_base

        if "basis" in siesta_inputs:
            new_basis = siesta_inputs["basis"].clone()
        else:
            new_basis = orm.Dict(dict={})

        new_basis["%block pao-basis"] = pao_block

        siesta_inputs["basis"] = new_basis

        run_extra = self.submit(SiestaBaseWorkChain, **siesta_inputs)

        self.report(
            f'Launched SiestaBaseWorkChain<{run_extra.pk}> with additional orbitals specified in `add_orbital`.'
        )

        return ToContext(extra_run=run_extra)

    def check_extra(self):
        """
        Check the outputs of the extra run with additional orbitals
        """
        extra_run = self.ctx.extra_run
        monitoring_quantity = self.ctx.monitoring_quantity

        if not extra_run.is_finished_ok:
            return self.exit_codes.ERROR_EXTRA_RUN

        self.report('Concluded calculations to test additional_orbitals. Analyzing results.')

        monitoring_quantity_value = extra_run.outputs.output_parameters[monitoring_quantity]

        if monitoring_quantity_value <= self.ctx.good_sizes_calc.outputs.output_parameters[monitoring_quantity]:
            self.ctx.good_sizes_calc = extra_run
            self.report(f'Additional orbitals led to a lower {monitoring_quantity}')
        else:
            self.report(f'Additional orbitals did not lead to a lower {monitoring_quantity}')

    def run_optimizer(self):
        """
        basis size. A check is than on what basis size gives the smallest "sizes_monitored_quantity".
        This quantity can be defined in input. By default it is the one used for the simplex and
        specified in the input "simplex.output_name".
        """

        if "ion_files" in self.ctx.good_sizes_calc.get_outgoing().nested():
            ions = self.ctx.good_sizes_calc.get_outgoing().nested()["ion_files"]

        global_enrgyshift = self.inputs.optimization_schema.global_energy_shift.value
        global_splitnorm = self.inputs.optimization_schema.global_split_norm.value
        charge_conf = self.inputs.optimization_schema.charge_confinement.value

        basis_dict, variables = extract_pao_block(global_splitnorm, global_enrgyshift, charge_conf, **ions)

        siesta_inputs = self.exposed_inputs(SimplexBasisOptimization).siesta_base

        if "basis" in siesta_inputs:
            new_basis = siesta_inputs["basis"].clone()
            translated_dict = FDFDict(new_basis.get_dict())
            for key in translated_dict:
                if key in ("paosplitnorm", "paoenergyshift"):
                    new_basis.attributes.pop(translated_dict.get_last_untranslated_key(key), None)
            new_basis.update_dict(basis_dict)
        else:
            new_basis = orm.Dict(dict=basis_dict)

        siesta_inputs["basis"] = new_basis

        simplex_inputs = self.exposed_inputs(SimplexBasisOptimization).simplex
        simplex_inputs["variables_dict"] = orm.Dict(dict=variables)
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
