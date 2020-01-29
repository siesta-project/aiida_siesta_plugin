# -*- coding: utf-8 -*-

from __future__ import absolute_import
import six
from six.moves import map

from aiida import orm
from aiida.engine import exceptions
from aiida.common import EntryPointError
from aiida.common.lang import override
from aiida.engine import CalcJob, WorkChain, ToContext, append_, while_
from aiida.plugins import CalculationFactory
from aiida.common import AttributeDict, AiidaException
from aiida.plugins.entry_point import get_entry_point_names, load_entry_point
from aiida_siesta.data.common import get_pseudos_from_structure
from aiida_siesta.calculations.siesta import SiestaCalculation

from .utils import register_error_handler, ErrorHandlerReport, ErrorHandler


class UnexpectedCalculationFailure(AiidaException):
    """Raised when a calculation has failed for an unknown reason."""
    pass


SiestaCalculation = CalculationFactory('siesta.siesta')


class SiestaBaseWorkChain(WorkChain):
    """
    Base Workchain to launch a total energy calculation via Siesta
    """

    _calculation_class = SiestaCalculation
    _error_handler_entry_point = 'aiida_siesta.workflow_error_handlers.base'

    def __init__(self, *args, **kwargs):
        super(SiestaBaseWorkChain, self).__init__(*args, **kwargs)

    #Next two functions are needed only if you specify an error handler registry and
    #you give it an entry point. It is not defined at the moment!
    @override
    def load_instance_state(self, saved_state, load_context):
        """Load the process instance from a saved state.
        :param saved_state: saved state of existing process instance
        :param load_context: context for loading instance state
        """
        super(SiestaBaseWorkChain,self).load_instance_state(saved_state, load_context)
        self._load_error_handlers()

    def _load_error_handlers(self):
        # If an error handler entry point is defined, load them. If the plugin cannot be loaded log it and pass
        if self._error_handler_entry_point is not None:
            for entry_point_name in get_entry_point_names(self._error_handler_entry_point):
                try:
                    load_entry_point(self._error_handler_entry_point, entry_point_name)
                    self.logger.info(
                        "loaded the '{}' entry point for the '{}' error handlers category"
                        .format(entry_point_name, self._error_handler_entry_point, plugin))
                except EntryPointError as exception:
                    self.logger.warning(
                        "failed to load the '{}' entry point for the '{}' error handlers: '{}'"
                        .format(entry_point_name, self._error_handler_entry_point, exception))


    @classmethod
    def define(cls, spec):
        super(SiestaBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=orm.Code)
        spec.input('structure', valid_type=orm.StructureData)
        spec.input_namespace('pseudos', required=False, dynamic=True)
        spec.input('pseudo_family', valid_type=orm.Str, required=False)
        spec.input('parent_folder', valid_type=orm.RemoteData, required=False)
        spec.input('kpoints', valid_type=orm.KpointsData, required=False)
        spec.input('bandskpoints', valid_type=orm.KpointsData, required=False)
        #Required by the plugin
        spec.input('parameters', valid_type=orm.Dict)
        spec.input('basis', valid_type=orm.Dict, required=False)
        spec.input('settings', valid_type=orm.Dict, required=False)
        #Required by any CalcJob
        spec.input('options', valid_type=orm.Dict)

        spec.input(
            'max_iterations',
            valid_type=orm.Int,
            default=orm.Int(5),
            help=
            'maximum number of iterations the workchain will restart the calculation to finish successfully'
        )
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
            default=orm.Bool(False),
            help=
            'if True, work directories of all called calculation will be cleaned at the end of execution'
        )

        spec.outline(
            cls.setup,
            cls.validate_pseudo_potentials,
            while_(cls.should_run_siesta)(
                cls.run_siesta,
                cls.inspect_siesta,
            ),
            cls.run_results,
        )

        # commented out, since it's more clear to provide list of outputs explicitly:
        # spec.dynamic_output()

        spec.output('forces_and_stress',
                    valid_type=orm.ArrayData,
                    required=False)
        spec.output('bands', valid_type=orm.BandsData, required=False)
        spec.output('output_structure',
                    valid_type=orm.StructureData,
                    required=False)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('remote_folder', valid_type=orm.RemoteData)

        spec.exit_code(
            100,
            'ERROR_ITERATION_RETURNED_NO_CALCULATION',
            message=
            'the run_calculation step did not successfully add a calculation node to the context'
        )
        spec.exit_code(101,
                       'ERROR_MAXIMUM_ITERATIONS_EXCEEDED',
                       message='the maximum number of iterations was exceeded')
        spec.exit_code(
            102,
            'ERROR_UNEXPECTED_CALCULATION_STATE',
            message=
            'the calculation finished with an unexpected calculation state')
        spec.exit_code(
            103,
            'ERROR_SECOND_CONSECUTIVE_SUBMISSION_FAILURE',
            message='the calculation failed to submit, twice in a row')
        spec.exit_code(
            104,
            'ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE',
            message=
            'the calculation failed for an unknown reason, twice in a row')

        spec.exit_code(
            301,
            'ERROR_INVALID_INPUT_PSEUDO_POTENTIALS',
            message=
            "the explicitly passed 'pseudos' or 'pseudo_family' could not be used to get the necessary potentials"
        )

        spec.exit_code(501,
                       'ERROR_WORKFLOW_FAILED',
                       message="Workflow did not succeed")

    def setup(self):
        """
        Initialize context variables
        """
        #self.report("Entering setup in Base Workchain")

        self.ctx.calc_name = 'SiestaCalculation'
        self.ctx.unexpected_failure = False
        self.ctx.submission_failure = False
        self.ctx.max_iterations = self.inputs.max_iterations.value
        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.ctx.iteration = 0
        #
        self.ctx.scf_did_not_converge = False
        self.ctx.geometry_did_not_converge = False
        self.ctx.want_band_structure = False
        self.ctx.out_of_time = False

        # Define inputs for the context of the workflow
        # the orm.Dict are back to be dictionaries because in this
        # way one can modify them. prepare_process_inputs reverts
        # them to orm.Dict
        # The few mandatory inputs (pseudo is handled later):
        self.ctx.inputs = {
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'pseudos': {},
            'parameters': self.inputs.parameters.get_dict(),
            'metadata': {
                'options': self.inputs.options.get_dict(),
            }
        }
        # Now the optional inputs
        if 'kpoints' in self.inputs:
            self.ctx.inputs['kpoints'] = self.inputs.kpoints
        if 'basis' in self.inputs:
            self.ctx.inputs['basis'] = self.inputs.basis.get_dict()
        if 'settings' in self.inputs:
            self.ctx.inputs['settings'] = self.inputs.settings.get_dict()
        if 'bandskpoints' in self.inputs:
            self.ctx.want_band_structure = True
            self.ctx.inputs['bandskpoints'] = self.inputs.bandskpoints

        #TO_DO soon, restart!
        # if 'parent_folder' in self.inputs:
        #     self.ctx.has_parent_folder = True
        #     self.ctx.inputs['parent_folder'] = self.inputs.parent_folder

        # Prevent SiestaCalculation from being terminated by scheduler
        max_wallclock_seconds = self.ctx.inputs['metadata']['options'][
            'max_wallclock_seconds']
        self.ctx.inputs['parameters']['max-walltime'] = max_wallclock_seconds

        return

    def prepare_pseudo_inputs(self,
                              structure,
                              pseudos=None,
                              pseudo_family=None):
        from aiida.orm import Str

        if pseudos and pseudo_family:
            raise ValueError(
                'you cannot specify both "pseudos" and "pseudo_family"')
        elif pseudos is None and pseudo_family is None:
            raise ValueError(
                'neither an explicit pseudos dictionary nor a pseudo_family was specified'
            )
        elif pseudo_family:
            # This will already raise some exceptions, potentially, like the ones below
            pseudos = get_pseudos_from_structure(structure,
                                                 pseudo_family.value)
        elif isinstance(pseudos, (six.string_types, Str)):
            raise TypeError(
                'you passed "pseudos" as a string - maybe you wanted to pass it as "pseudo_family" instead?'
            )

        for kind in structure.get_kind_names():
            if kind not in pseudos:
                raise ValueError(
                    'no pseudo available for element {}'.format(kind))

        return pseudos

    def validate_pseudo_potentials(self):
        """
        Validate the inputs related to pseudopotentials to check that we have the minimum required
        amount of information to be able to run a SiestaCalculation
        """
        structure = self.inputs.structure
        pseudos = self.inputs.get('pseudos', None)
        pseudo_family = self.inputs.get('pseudo_family', None)

        try:
            self.ctx.inputs['pseudos'] = self.prepare_pseudo_inputs(
                structure, pseudos, pseudo_family)
        except ValueError as exception:
            self.report('{}'.format(exception))
            return self.exit_codes.ERROR_INVALID_INPUT_PSEUDO_POTENTIALS

    def should_run_siesta(self):
        """
        Return whether a siesta restart calculation should be run, which
        is the case as long as the last calculation was not converged
        successfully and the maximum number of restarts has not yet
        been exceeded
        """
        return ((not self.ctx.is_finished)
                and (self.ctx.iteration < self.ctx.max_iterations))

    def prepare_process_inputs(self, inputs):
        """
        Prepare the inputs for submission for the given process, according
        to its spec. That is to say that when an input is found in the
        inputs that corresponds to an input port in the spec of the
        process that expects a Dict, yet the value in the inputs is a
        plain dictionary, the value will be wrapped in by the Dict
        class to create a valid input.

        :param process: sub class of Process for which to prepare the inputs dictionary
        :param inputs: a dictionary of inputs intended for submission of the process
        :return: a dictionary with all bare dictionaries wrapped in Dict if dictated by process spec

        """
        prepared_inputs = AttributeDict()
        process_spec = SiestaCalculation.spec()

        for key, value in six.iteritems(inputs):

            if key not in process_spec.inputs:
                continue

            if process_spec.inputs[key].valid_type == orm.Dict and isinstance(
                    value, dict):
                prepared_inputs[key] = orm.Dict(dict=value)
            else:
                prepared_inputs[key] = value

        return prepared_inputs

    def run_siesta(self):
        """
        Run a new SiestaCalculation or restart from a previous
        SiestaCalculation run in this workchain

        """

        self.ctx.iteration += 1

        # wrapping inputs to Dict if they are dicts, or returning raw
        try:
            wrapped_inputs = self.ctx.inputs
        except AttributeError:
            raise ValueError(
                'no calculation input dictionary was defined in self.ctx.inputs'
            )

        inputs = self.prepare_process_inputs(wrapped_inputs)
        calculation = self.submit(SiestaCalculation, **inputs)
        self.report('launching {}<{}> iteration #{}'.format(
            self.ctx.calc_name, calculation.pk, self.ctx.iteration))

        return ToContext(calculations=append_(calculation))

    def inspect_siesta(self):
        """
        Analyse the results of the previous SiestaCalculation, checking
        whether it finished successfully, or if not troubleshoot the
        cause and adapt the input parameters accordingly before
        restarting, or abort if unrecoverable error was found
        """
        try:
            calculation = self.ctx.calculations[self.ctx.iteration - 1]
        except IndexError:
            self.report('iteration {} finished without returning a {}'.format(
                self.ctx.iteration, self.ctx.calc_name))
            return self.exit_codes.ERROR_ITERATION_RETURNED_NO_CALCULATION

        exit_code = None

        # Done: successful completion of last calculation
        if calculation.is_finished_ok:
            #self.report('{}<{}> completed successfully'
            #            .format(self.ctx.calc_name, calculation.pk))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True

        # Abort: exceeded maximum number of retries
        elif self.ctx.iteration >= self.inputs.max_iterations.value:
            self.report(
                'reached the maximumm number of iterations {}: last ran {}<{}>'
                .format(self.inputs.max_iterations.value, self.ctx.calc_name,
                        calculation.pk))
            exit_code = self.exit_codes.ERROR_MAXIMUM_ITERATIONS_EXCEEDED

        # Retry or abort: calculation finished or failed
        else:

            # Calculation was at least submitted successfully, so we reset the flag
            self.ctx.submission_failure = False

            # calculation failed, try to salvage it or handle any unexpected failures
            try:
                exit_code = self._handle_calculation_failure(calculation)
            except UnexpectedCalculationFailure as exception:
                exit_code = self._handle_unexpected_failure(
                    calculation, exception)
                self.ctx.unexpected_failure = True

        return exit_code

    def run_results(self):
        """
        Attach the output parameters and retrieved folder of the last
        calculation to the outputs

        """


        for name, port in six.iteritems(self.spec().outputs):

            try:
                node = self.ctx.restart_calc.get_outgoing(
                    link_label_filter=name).one().node
            except ValueError:
                if port.required:
                    self.report(
                        "the process spec specifies the output '{}' as required but was not an output of {}<{}>"
                        .format(name, self.ctx.calc_name,
                                self.ctx.restart_calc.pk))
            else:
                self.out(name, node)
                #self.report("attaching the node {}<{}> as '{}'"
                #            .format(node.__class__.__name__, node.pk, name))

        self.report('Base workchain completed after {} iterations'.format(
            self.ctx.iteration))
    def on_terminated(self):
        """
        If the clean_workdir input was set to True, recursively collect
        all called Calculations by ourselves and our called
        descendants, and clean the remote folder for the CalcJobNode
        instances

        """
        super(SiestaBaseWorkChain, self).on_terminated()

        if self.inputs.clean_workdir.value is False:
            return

        cleaned_calcs = []

        for called_descendant in self.calc.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                try:
                    called_descendant.outputs.remote_folder._clean()
                    cleaned_calcs.append(called_descendant.pk)
                except (IOError, OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report('cleaned remote folders of calculations: {}'.format(
                ' '.join(map(str, cleaned_calcs))))

    def _handle_submission_failure(self, calculation):
        """
        The submission of the calculation has failed. If the
        submission_failure flag is set to true, this is the second
        consecutive submission failure and we abort the workchain
        Otherwise we restart once more.

        """
        if self.ctx.submission_failure:
            self.report(
                'submission for {}<{}> failed for the second consecutive time'.
                format(self.ctx.calc_name, calculation.pk))
            return self.exit_codes.ERROR_SECOND_CONSECUTIVE_SUBMISSION_FAILURE

        else:
            self.report(
                'submission for {}<{}> failed, restarting once more'.format(
                    self.ctx.calc_name, calculation.pk))

    def _handle_unexpected_failure(self, calculation, exception=None):
        """
        The calculation has failed for an unknown reason and could not be
        handled. If the unexpected_failure flag is true, this is the
        second consecutive unexpected failure and we abort the
        workchain.  Otherwise we restart once more.

        """
        if exception:
            self.report('{}'.format(exception))

        if self.ctx.unexpected_failure:
            self.report(
                'failure of {}<{}> could not be handled for a second consecutive time'
                .format(self.ctx.calc_name, calculation.pk))
            return self.exit_codes.ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE

        else:
            self.report(
                'failure of {}<{}> could not be handled, restarting once more'.
                format(self.ctx.calc_name, calculation.pk))

    def _handle_calculation_failure(self, calculation):
        """
        The calculation has failed so we try to analyze the reason and
        change the inputs accordingly for the next calculation. If the
        calculation failed, but did so cleanly, we set it as the
        restart_calc, in all other cases we do not replace the
        restart_calc

        """
        try:
            outputs = calculation.outputs.output_parameters.get_dict(
            )['warnings']
            # _ = outputs['warnings']
            # _ = outputs['parser_warnings']
        except (AttributeError, KeyError) as exception:
            raise UnexpectedCalculationFailure(exception)

        is_handled = False
        handler_report = None

        # Sort the handlers based on their priority in reverse order
        handlers = sorted(self._error_handlers,
                          key=lambda x: x.priority,
                          reverse=True)

        if not handlers:
            raise UnexpectedCalculationFailure(
                'no calculation error handlers were registered')

        for handler in handlers:
            # print(handler)
            handler_report = handler.method(self, calculation)

            # If at least one error is handled, we consider the
            # calculation failure handled
            if handler_report and handler_report.is_handled:
                is_handled = True

            # After certain error handlers, we may want to skip all
            # other error handling
            if handler_report and handler_report.do_break:
                break

        # If none of the executed error handlers reported that they
        # handled an error, the failure reason is unknown
        if not is_handled:
            raise UnexpectedCalculationFailure(
                'calculation failure was not handled')

        # The last called error handler may not necessarily have
        # returned a handler report
        if handler_report:
            return handler_report.exit_code

        return


@register_error_handler(SiestaBaseWorkChain, 130)
def _handle_error_geom_not_conv(self, calculation):
    """
    At the end of the scf cycle, the geometry convergence was not
    reached.  We need to restart from the previous calculation
    """

    self.report(
        'SiestaCalculation<{}> did not reach geometry convergence. Will restart.'
        .format(calculation.pk))

    g = calculation
    # We need to take care here of passing the
    # output geometry of old_calc to the new calculation
    if g.outputs.output_parameters.attributes["variable_geometry"]:
        self.ctx.inputs['structure'] = g.outputs.output_structure

    #The most important line. The presence of
    #parent_calc_folder triggers the real restart
    #meaning the copy of the .DM and the
    #addition of use-saved-dm to the parameters

    self.ctx.inputs['parent_calc_folder'] = g.outputs.remote_folder

    self.ctx.restart_calc = calculation

    return ErrorHandlerReport(True, False)


@register_error_handler(SiestaBaseWorkChain, 120)
def _handle_error_scf_not_conv(self, calculation):
    """
    SCF convergence was not reached.  We need to restart from the
    previous calculation without changing any of the input parameters.
    """

    self.report(
        'SiestaCalculation<{}> did not achieve scf convergence. Will restart.'.
        format(calculation.pk))

    # The most important line. The presence of
    # parent_calc_folder triggers the real restart
    # meaning the copy of the .DM and the
    # addition of use-saved-dm to the parameters

    self.ctx.inputs['parent_calc_folder'] = calculation.outputs.remote_folder

    self.ctx.restart_calc = calculation

    return ErrorHandlerReport(True, False)
