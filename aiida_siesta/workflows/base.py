# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Bool, Int, Str
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, while_

from aiida_siesta.data.psf import PsfData, get_pseudos_from_structure
from aiida_siesta.calculations.siesta import SiestaCalculation


class SiestaBaseWorkChain(WorkChain):
    """
    Base Workchain to launch a total energy calculation via Siesta
    """
    def __init__(self, *args, **kwargs):
        super(SiestaBaseWorkChain, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(SiestaBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input_group('pseudos', required=False)
        spec.input('pseudo_family', valid_type=Str, required=False)
        spec.input('parent_folder', valid_type=RemoteData, required=False)
        spec.input('kpoints', valid_type=KpointsData)
        spec.input('bandskpoints', valid_type=KpointsData, required=False)
        spec.input('parameters', valid_type=ParameterData)
        spec.input('basis', valid_type=ParameterData)
        spec.input('settings', valid_type=ParameterData)
        spec.input('options', valid_type=ParameterData)
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
        spec.input('max_iterations', valid_type=Int, default=Int(10))
        spec.outline(
            cls.setup,
            cls.validate_pseudo_potentials,
            while_(cls.should_run_siesta)(
                cls.run_siesta,
                cls.inspect_siesta,
            ),
            cls.run_results,
        )
        spec.dynamic_output()

    def setup(self):
        """
        Initialize context variables
        """
        self.report("Entering setup in Base Workchain")
        
        self.ctx.max_iterations = self.inputs.max_iterations.value
        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.ctx.iteration = 0
        #
        self.ctx.scf_did_not_converge = False
        self.ctx.geometry_did_not_converge = False
        self.ctx.want_band_structure = False
        self.ctx.out_of_time = False

        # Define convenience dictionary of inputs for SiestaCalculation
        self.ctx.inputs = {
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'pseudo': {},
            'kpoints': self.inputs.kpoints,
            'parameters': self.inputs.parameters.get_dict(),
            'basis': self.inputs.basis.get_dict(),
            'settings': self.inputs.settings.get_dict(),
            '_options': self.inputs.options.get_dict(),
        }
        if 'bandskpoints' in self.inputs:
            self.ctx.want_band_structure = True
            self.ctx.inputs['bandskpoints'] = self.inputs.bandskpoints
        #
        #  This is of limited usefulness, as would need
        #  to check that the items in the folder are compatible
        #  with the rest of the calculation's parameters
        
        if 'parent_folder' in self.inputs:
            self.ctx.has_parent_folder = True
            self.ctx.inputs['parent_folder'] = self.inputs.parent_folder

        # Prevent SiestaCalculation from being terminated by scheduler
        max_wallclock_seconds = self.ctx.inputs['_options']['max_wallclock_seconds']
        self.ctx.inputs['parameters']['max-walltime'] = max_wallclock_seconds

        return

    def validate_pseudo_potentials(self):
        """
        Validate the inputs related to pseudopotentials to check that we have the minimum required
        amount of information to be able to run a SiestaCalculation
        """

        if all([key not in self.inputs for key in ['pseudos', 'pseudo_family']]):
            self.abort_nowait('neither explicit pseudos nor a pseudo_family was specified in the inputs')
            return
        elif all([key in self.inputs for key in ['pseudos', 'pseudo_family']]):
            self.report('both explicit pseudos as well as a pseudo_family were specified: using explicit pseudos')
            self.ctx.inputs['pseudo'] = self.inputs.pseudos
        elif 'pseudos' in self.inputs:
            self.report('only explicit pseudos were specified: using explicit pseudos')
            self.ctx.inputs['pseudo'] = self.inputs.pseudos
        elif 'pseudo_family' in self.inputs:
            self.report('only a pseudo_family was specified: using pseudos from pseudo_family')
            structure = self.inputs.structure
            self.ctx.inputs['pseudo'] = get_pseudos_from_structure(structure, self.inputs.pseudo_family.value)

        for kind in self.inputs.structure.get_kind_names():
            if kind not in self.ctx.inputs['pseudo']:
                self.abort_nowait('no pseudo available for element {}'.format(kind))
            elif not isinstance(self.ctx.inputs['pseudo'][kind], PsfData):
                self.abort_nowait('pseudo for element {} is not of type PsfData'.format(kind))

    def should_run_siesta(self):
        """
        Return whether a siesta restart calculation should be run, which
        is the case as long as the last calculation was not converged
        successfully and the maximum number of restarts has not yet
        been exceeded
        """
        return ( (not self.ctx.is_finished) and (self.ctx.iteration < self.ctx.max_iterations) )

    def run_siesta(self):
        """
        Run a new SiestaCalculation or restart from a previous
        SiestaCalculation run in this workchain

        """
        self.report("Running Siesta")
        
        self.ctx.iteration += 1

        # Create local copy of general initial inputs stored in the context
        # and adapt for next calculation
        
        local_inputs = dict(self.ctx.inputs)

        # Indicates if restarting or doing calculations from scratch.
        # left from the original WorkChain template for QE pw.x
        # TODO Clarify on that matter
        #
        # if self.ctx.iteration == 1 and 'parent_folder' in self.inputs:
        #  For Siesta, we should decide on which files to actually
        #  use in the parent_folder. It might be enough to specify
        #  'dm-use-save-DM'
        #
        #     inputs['parameters']['CONTROL']['restart_mode'] = 'restart'
        #     inputs['parent_folder'] = self.inputs.parent_folder
        # elif self.ctx.restart_calc:
        #     inputs['parameters']['CONTROL']['restart_mode'] = 'restart'
        #     inputs['parent_folder'] = self.ctx.restart_calc.out.remote_folder
        # else:
        #     inputs['parameters']['CONTROL']['restart_mode'] = 'from_scratch'

        # NOTE really the logic should be here

        if self.ctx.restart_calc:
            local_inputs['parent_folder'] = self.ctx.restart_calc.out.remote_folder

        if self.ctx.scf_did_not_converge:
            local_inputs['parameters']['dm-use-save-dm'] = True
            self.report('Re-using previous DM')

        # We need to add here the previous structure, for cases of
        # geometry optimization
        #
        if self.ctx.geometry_did_not_converge:
            
            # The previous calculation, even if failed, should
            # have produced an 'output_structure' node containing the
            # last recorded geometry in the XML file.
            # 
            # So we benefit from all the geometry optimization work
            # done in that calculation (except from the internal state
            # of the optimizer...)
            #
            # Another route to geometry re-use is to employ the
            # information in files that checkpoint the geometry:
            #
            # 'move XXX.STRUCT_OUT TO XXX.STRUCT_IN'
            #  and add the 'use-struct-in' fdf option
            # 
            local_inputs['structure'] = self.ctx.restart_calc.out.output_structure
            self.report('Re-using previous output_structure')
            # --- maybe decide whether to actually use the DM... or to extrapolate...
            local_inputs['parameters']['dm-use-save-dm'] = True
            self.report('Re-using previous DM')

        
        local_inputs['parameters'] = ParameterData(dict=local_inputs['parameters'])

        local_inputs['basis'] = ParameterData(dict=local_inputs['basis'])
        local_inputs['settings'] = ParameterData(dict=local_inputs['settings'])

        process = SiestaCalculation.process()
        running = submit(process, **local_inputs)

        self.report('launching SiestaCalculation<{}> iteration #{}'.format(running.pid, self.ctx.iteration))

        return ToContext(calculation=running)

    def inspect_siesta(self):

        """
        Analyse the results of the previous SiestaCalculation, checking
        whether it finished successfully, or if not troubleshoot the
        cause and adapt the input parameters accordingly before
        restarting, or abort if unrecoverable error was found
        """
        try:
            calculation = self.ctx.calculation
        except Exception:
            self.abort_nowait('The previous iteration finished without returning a SiestaCalculation')
            return

        expected_states = [calc_states.FINISHED, calc_states.FAILED, calc_states.SUBMISSIONFAILED]

        # Done: successful convergence of last calculation
        if calculation.has_finished_ok():
            self.report('converged successfully after {} iterations'.format(self.ctx.iteration))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True

        # Abort: exceeded maximum number of retries
        elif self.ctx.iteration >= self.ctx.max_iterations:
            self.report('reached the max number of iterations {}'.format(self.ctx.max_iterations))
            self.abort_nowait('last ran SiestaCalculation<{}>'.format(calculation.pk))

        # Abort: unexpected state of last calculation
        elif calculation.get_state() not in expected_states:
            self.abort_nowait('unexpected state ({}) of SiestaCalculation<{}>'.format(
                calculation.get_state(), calculation.pk))

        # Retry: submission failed, try to restart or abort
        # NOTE This handler is not implemented
        elif calculation.get_state() in [calc_states.SUBMISSIONFAILED]:
            self._handle_submission_failure(calculation)

        # Retry: calculation failed
        # The FAILED state is the one actually reported for
        # non-converged calculations when the
        # 'scf-must-converge'
        # and/or
        # 'geometry-must-converge'
        # fdf options are specified. There are then 'FATAL'
        # lines in the MESSAGES file (and thus in the 'warnings'
        # list)
        #
        # We might also check 'FINISHED' calculations for
        # 'WARNING' SCF_NOT_CONV and/or GEOM_NOT_CONV lines.
        #
        elif calculation.get_state() in [calc_states.FAILED]:
            self._handle_calculation_failure(calculation)

        else:
            self.abort_nowait('This place should not be reached...')
        return

    def run_results(self):
        """
        Attach the output parameters and retrieved folder of the last calculation to the outputs
        """
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))
        self.out('output_parameters', self.ctx.restart_calc.out.output_parameters)
        self.out('remote_folder', self.ctx.restart_calc.out.remote_folder)
        self.out('retrieved', self.ctx.restart_calc.out.retrieved)

        if 'output_structure' in self.ctx.restart_calc.out:
            self.out('output_structure', self.ctx.restart_calc.out.output_structure)
        if 'bands_array' in self.ctx.restart_calc.out:
            self.out('bands_array', self.ctx.restart_calc.out.bands_array)

    def _handle_submission_failure(self, calculation):

        """
        The submission of the calculation has failed, if it was the second
        consecutive failure we abort the workchain, else we set the
        has_submission_failed flag and try again
        """
        self.abort_nowait('submission failed for the {} in iteration {}, but error handling is not implemented yet'
            .format(SiestaCalculation.__name__, self.ctx.iteration))

    def _handle_calculation_failure(self, calculation):
        """
        The calculation has failed so we try to analyze the reason and
        change the inputs accordingly for the next calculation. If the
        calculation failed, but did so cleanly, we set it as the
        restart_calc, in all other cases we do not replace the
        restart_calc

        """
        # Typical contents of the warnings list for a FAILED calculation:
        #
        #        "warnings": [
        #        "FATAL: GEOM_NOT_CONV: Geometry relaxation not converged",
        #        "FATAL: ABNORMAL_TERMINATION"
        
        warnings_list = calculation.out.output_parameters.get_dict()['warnings']

        # Formally we should be checking also for an OUT_OF_TIME fatal message,
        # but in this case Siesta attaches SCF or GEOM 'WARNING' lines, so the
        # checks below will cover it:
        #
        #  WARNING: SCF_NOT_CONV: SCF did not converge at wall time exhaustion
        #   (info): Geom step, scf iteration, dmax:   4  4    0.000413
        #  FATAL: OUT_OF_TIME: Time is up.
        
        # We should, however, report it:
        
        self.ctx.out_of_time = False
        for line in warnings_list:
            if u'FATAL: OUT_OF_TIME' in line:
                self.ctx.out_of_time = True
                self.report('Out of time in SiestaCalculation<{}>'.format(calculation.pk))

        # Note again that we check for the strings themselves, and not
        # for 'FATAL' or 'WARNING' qualifiers
        
        self.ctx.geometry_did_not_converge = False
        for line in warnings_list:
            if u'GEOM_NOT_CONV' in line:
                self.ctx.geometry_did_not_converge = True

        self.ctx.scf_did_not_converge = False
        for line in warnings_list:
            if u'SCF_NOT_CONV' in line:
                self.ctx.scf_did_not_converge = True

        # We might have run out of time during the analysis stage, which
        # includes the bands calculation
        
        self.ctx.restart_calc = calculation
                

    # def on_stop(self):
    #     """Clean remote folders of the SiestaCalculations that were run if
    #     the clean_workdir parameter was set to true in the Workchain
    #     inputs

    #     """
    #     super(SiestaBaseWorkChain, self).on_stop()

    #     if not self.inputs.clean_workdir.value:
    #         self.report('SiestaBase: remote folders will not be cleaned')
    #         return

    #     for calc in self.ctx.calculation:
    #         try:
    #             calc.out.remote_folder._clean()
    #             self.report('cleaned remote folder of {}<{}>'.format(calc.__class__.__name__, calc.pk))
    #         except Exception:
    #             pass
