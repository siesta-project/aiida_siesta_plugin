# -*- coding: utf-8 -*-
"""
Base workchain to run siesta calculation.

It handles error and it is based on the RestartWorkChain of aiida core.
"""

from aiida import orm
from aiida.common import AttributeDict
from aiida.common.exceptions import NotExistent
from aiida.engine import BaseRestartWorkChain, ProcessHandlerReport, process_handler, while_

from aiida_siesta.calculations.siesta import SiestaCalculation, bandskpoints_warnings, internal_structure
from aiida_siesta.utils.tkdict import FDFDict


def validate_options(value, _):
    """
    Validate options input port.
    """
    if value:
        if "max_wallclock_seconds" not in value.get_dict():
            return "the `max_wallclock_seconds` key is required in the options dict."


def validate_ps_fam(value, _):
    """
    Validate pseudo_family input port.
    """
    if value:
        try:
            orm.Group.get(label=value.value)
        except NotExistent:
            return f"{value.value} does not correspond to any known pseudo family."


def validate_inputs(value, _):
    """
    Validate the entire input namespace.

    It takes care to ckeck the consistency and compatibility of the inputed basis, pseudos, pseudofamilies and ions.
    Also calls the `bandskpoints_warnings` that issues warning about bandskpoints selection.
    It is similar to the `validate_inputs` of SiestaCalculation but with the additional complexity
    due to the presence of 'pseudo_family' input.
    """
    import warnings

    bandskpoints_warnings(value)

    if 'structure' in value:
        if 'basis' in value:
            structure = internal_structure(value["structure"], value["basis"].get_dict())
            if structure is None:
                return "Not possibe to specify `floating_sites` (ghosts) with the same name of a structure kind."
        else:
            structure = value["structure"]

        #Check each kind in the structure (including freshly added ghosts) have a corresponding pseudo or ion
        kinds = [kind.name for kind in structure.kinds]
        if 'ions' in value:
            quantity = 'ions'
            if 'pseudos' in value or 'pseudo_family' in value:
                warnings.warn("At least one ion file in input, all the pseudos or pseudo_family will be ignored")
        else:
            quantity = 'pseudos'
            if 'pseudos' not in value and 'pseudo_family' not in value:
                return "No `pseudos`, nor `ions`, nor `pseudo_family` specified in input"
            if 'pseudos' in value and 'pseudo_family' in value:
                return "You cannot specify both `pseudos` and `pseudo_family`"

        if 'pseudo_family' in value:
            group = orm.Group.get(label=value['pseudo_family'].value)
            try:
                group.get_pseudos(structure=structure)
            except ValueError:
                return "The pseudo family does not incude all the required pseudos"
        else:
            if set(kinds) != set(value[quantity].keys()):
                ps_io = ', '.join(list(value[quantity].keys()))
                kin = ', '.join(list(kinds))
                string_out = (
                    'mismatch between defined pseudos/ions and the list of kinds of the structure\n' +
                    f' pseudos/ions: {ps_io} \n kinds(including ghosts): {kin}'
                )
                return string_out


class SiestaBaseWorkChain(BaseRestartWorkChain):
    """
    Base Workchain to launch a total energy calculation via Siesta.
    """

    _process_class = SiestaCalculation
    _proc_exit_cod = _process_class.exit_codes

    @classmethod
    def define(cls, spec):
        """
        Define the specs.
        """
        super().define(spec)
        spec.expose_inputs(SiestaCalculation, exclude=('metadata',))
        spec.input('pseudo_family', valid_type=orm.Str, required=False, validator=validate_ps_fam)
        spec.input('options', valid_type=orm.Dict, validator=validate_options)

        spec.outline(
            cls.preprocess,
            cls.setup,
            cls.prepare_inputs,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
            cls.postprocess,
        )

        spec.inputs.validator = validate_inputs

        spec.expose_outputs(SiestaCalculation)

        spec.exit_code(403, 'ERROR_BASIS_POL', message='Basis polarization problem.')
        spec.exit_code(404, 'ERROR_BANDS_PARSING', message='Error in the parsing of bands')

    def preprocess(self):
        """
        Here a higher level WorkChain could put preprocesses.
        """

    def prepare_inputs(self):
        """
        Initialize context variables.

        Note that in context we must include exactly the same data nodes obtained in input (except the parameters
        where we add `max-walltime`)
        On the contrary, useless nodes would be created.
        """
        self.report("Preparing inputs of the SiestaBaseWorkChain")

        self.ctx.inputs = AttributeDict(self.exposed_inputs(SiestaCalculation))
        self.ctx.inputs['metadata'] = {'options': self.inputs.options.get_dict()}

        structure = self.inputs.structure

        if "pseudo_family" in self.inputs:
            fam_name = self.inputs.pseudo_family.value
            group = orm.Group.get(label=fam_name)
            if 'basis' in self.inputs:
                temp_structure = internal_structure(structure, self.inputs.basis.get_dict())
                self.ctx.inputs['pseudos'] = group.get_pseudos(structure=temp_structure)
            else:
                self.ctx.inputs['pseudos'] = group.get_pseudos(structure=structure)

    def postprocess(self):
        """
        Attach the output_namespaces to outputs.

        In theory, the BaseRestartWorkChain should already return all the output
        requested in spec if they are returned (output nodes) by the last completed
        process. However this fails for `output_namespaces` (issue #4623 aiida-core).
        For this reason I do the procedure to attach the `output_namespaces` here.
        """
        if "ion_files" in self.spec().outputs:
            node = self.ctx.children[self.ctx.iteration - 1]
            if "ion_files" in node.get_outgoing().nested():
                self.out("ion_files", node.get_outgoing().nested()["ion_files"])

    @process_handler(priority=70, exit_codes=_proc_exit_cod.GEOM_NOT_CONV)  #pylint: disable = no-member
    def handle_error_geom_not_conv(self, node):
        """
        Geometry convergence was not reached.We need to restart from the previous calculation.
        """

        self.report(f'SiestaCalculation<{node.pk}> did not reach geometry convergence')

        # We need to take care of passing the output geometry of old_calc to the new calculation.
        if node.outputs.output_parameters.attributes["variable_geometry"]:
            self.ctx.inputs['structure'] = node.outputs.output_structure

        # The presence of `parent_calc_folder` triggers the real restart, so we add it.
        self.ctx.inputs['parent_calc_folder'] = node.outputs.remote_folder

        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=80, exit_codes=_proc_exit_cod.SCF_NOT_CONV)  #pylint: disable = no-member
    def handle_error_scf_not_conv(self, node):
        """
        SCF convergence was not reached.

        We need to restart from the previous calculation without changing any of the input parameters.
        """

        self.report(f'SiestaCalculation<{node.pk}> did not achieve scf convergence.')

        # We need to take care of passing the output geometry of old_calc to the new calculation.
        if node.outputs.output_parameters.attributes["variable_geometry"]:
            self.ctx.inputs['structure'] = node.outputs.output_structure

        # The presence of `parent_calc_folder` triggers the real restart, so we add it.
        self.ctx.inputs['parent_calc_folder'] = node.outputs.remote_folder

        #Should be also increase the number of scf max iterations?

        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=90, exit_codes=_proc_exit_cod.SPLIT_NORM)  #pylint: disable = no-member
    def handle_error_split_norm(self, node):
        """
        The split_norm parameter was too small.

        We need to change it and restart. The minimum split_norm is stored in the logs of the old calculation.
        This error happens only at the beginning of the run, therefore no real restart needed.
        Just a new calculation with a new split_norm.
        """

        self.report(f'SiestaCalculation<{node.pk}> crashed with split_norm issue.')

        #Retrive the minimum split norm from the logs of failed calc.
        logs = orm.Log.objects.get_logs_for(node)
        for log in logs:
            if "Error in split_norm option" in log.message:
                mylog = log.message.split()
        new_split_norm = float(mylog[-1]) + 0.001

        # We want to understand the presence of "pao-split-norm" in input and:
        # 1) if present, we change its value to the minimum allowed
        # 2) if not present, we activate pao-SplitTailNorm
        # As we don't know in which sintax the user passed "pao-split-norm (remember that every fdf variant
        # is allowed), we translate the original dict to a FDFDict that is aware of equivalent keyword.
        transl_basis = FDFDict(self.ctx.inputs["basis"].get_dict())
        glob_split_norm = False
        for key in transl_basis:
            if key == "paosplitnorm":
                glob_split_norm = True

        if glob_split_norm:
            self.report('Resetting the pao-split-norm global value')
            transl_basis["pao-split-norm"] = new_split_norm
        else:
            self.report('Adding pao-SplitTailNorm to solve the split_norm problem')
            transl_basis["pao-SplitTailNorm"] = True

        new_basis = orm.Dict(dict=transl_basis)
        self.ctx.inputs["basis"] = new_basis

        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=89, exit_codes=_proc_exit_cod.BASIS_POLARIZ)  #pylint: disable = no-member
    def handle_error_basis_pol(self, node):  #pylint: disable = unused-argument
        """
        Problem in generating plarized orbitals.

        For the moment, we don't handle this error, but we terminate the WorkChain with
        a specific error code.
        """
        return ProcessHandlerReport(True, self.exit_codes.ERROR_BASIS_POL)

    #pylint: disable = no-member
    @process_handler(priority=60, exit_codes=[_proc_exit_cod.BANDS_PARSE_FAIL, _proc_exit_cod.BANDS_FILE_NOT_PRODUCED])
    def handle_error_bands(self, node):  #pylint: disable = unused-argument
        """
        Error on bands.

        If an error in the parsing of bands occours in the SiestaCalculation (node here), we expose all
        the output ports node that have been produced (SiestaCalculation is designed to produce the
        output_parameter and stress/forcess port before the check on the bands outputs) and then we
        stop the workchain with a specific error code. We exclude the "retrieved" output port as
        it refers only to the underline calculation, not to the WorkChain itself.
        """
        for name in node.outputs:
            if name != "retrieved":
                output = node.get_outgoing(link_label_filter=name).one().node
                self.out(name, output)

        return ProcessHandlerReport(True, self.exit_codes.ERROR_BANDS_PARSING)

    @classmethod
    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
        """
        Get the input generator.
        """
        from aiida_siesta.utils.protocols_system.input_generators import BaseWorkChainInputGenerator
        return BaseWorkChainInputGenerator(cls)
