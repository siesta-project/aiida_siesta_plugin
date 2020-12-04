from aiida import orm
from aiida.engine import BaseRestartWorkChain, ProcessHandlerReport, process_handler, while_
from aiida_siesta.data.common import get_pseudos_from_structure
from aiida_siesta.calculations.siesta import SiestaCalculation


def prepare_pseudo_inputs(structure, pseudos, pseudo_family):

    if pseudos is not None and pseudo_family is not None:
        raise ValueError('you cannot specify both "pseudos" and "pseudo_family"')
    elif pseudos is None and pseudo_family is None:
        raise ValueError('neither an explicit pseudos dictionary nor a pseudo_family was specified')
    elif pseudo_family is not None:
        # This will already raise some exceptions
        pseudos = get_pseudos_from_structure(structure, pseudo_family.value)

    for kind in structure.get_kind_names():
        if kind not in pseudos:
            raise ValueError('no pseudo available for element {}'.format(kind))

    return pseudos


class SiestaBaseWorkChain(BaseRestartWorkChain):
    """
    Base Workchain to launch a total energy calculation via Siesta
    """
    _process_class = SiestaCalculation
    _proc_exit_cod = _process_class.exit_codes

    @classmethod
    def define(cls, spec):
        super(SiestaBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=orm.Code)
        spec.input('structure', valid_type=orm.StructureData)
        spec.input_namespace('pseudos', required=False, dynamic=True)
        spec.input('pseudo_family', valid_type=orm.Str, required=False)
        spec.input('parent_calc_folder', valid_type=orm.RemoteData, required=False)
        spec.input('kpoints', valid_type=orm.KpointsData, required=False)
        spec.input('bandskpoints', valid_type=orm.KpointsData, required=False)
        spec.input('parameters', valid_type=orm.Dict)
        spec.input('basis', valid_type=orm.Dict, required=False)
        spec.input('settings', valid_type=orm.Dict, required=False)
        spec.input('options', valid_type=orm.Dict)

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

        spec.output('forces_and_stress', valid_type=orm.ArrayData, required=False)
        spec.output('bands', valid_type=orm.BandsData, required=False)
        spec.output('output_structure', valid_type=orm.StructureData, required=False)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('remote_folder', valid_type=orm.RemoteData)

        spec.exit_code(403, 'ERROR_BASIS_POL', message='Basis polarization problem.')
        spec.exit_code(404, 'ERROR_BANDS_PARSING', message='Error in the parsing of bands')

    def preprocess(self):
        """
        Here a higher level WorkChain could put preprocesses
        """

    def prepare_inputs(self):
        """
        Initialize context variables
        """
        self.report("Preparing inputs of the SiestaBaseWorkChain")

        structure = self.inputs.structure
        pseudo_family = self.inputs.get('pseudo_family', None)
        #The port 'pseudos' is an 'input_namespace', therefore is never undefined in the
        #current aiida implementation (see Issue #142 plumpy), but it is an empty dictionary
        #if pseudos are not passed in input.
        #Therefore 'pseudos = self.inputs.get('pseudos', None)' never gives None. Better:
        pseudos = None
        if "pseudos" in self.inputs:  #in case in the future Issue #142 will be solved
            if self.inputs.pseudos:
                pseudos = self.inputs.pseudos

        self.ctx.inputs = {
            'code': self.inputs.code,
            'structure': structure,
            'pseudos': prepare_pseudo_inputs(structure, pseudos, pseudo_family),
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
        if 'parent_calc_folder' in self.inputs:
            self.ctx.inputs['parent_calc_folder'] = self.inputs.parent_calc_folder

        # Prevent SiestaCalculation from being terminated by scheduler
        max_wallclock_seconds = self.ctx.inputs['metadata']['options']['max_wallclock_seconds']
        self.ctx.inputs['parameters']['max-walltime'] = max_wallclock_seconds

        #Note: To pass pure dictionaries or orm.Dict is the same as the WorkChain
        #will take care of wrapping in orm.Dict the pure python dict before submission,
        #however this influences the way you fix problems in the hadlers above.

    def postprocess(self):
        """
        Here a higher level WorkChain could put postprocesses
        """

    @process_handler(priority=70, exit_codes=_proc_exit_cod.GEOM_NOT_CONV)  #pylint: disable = no-member
    def handle_error_geom_not_conv(self, node):
        """
        At the end of the scf cycle, the geometry convergence was not
        reached.  We need to restart from the previous calculation
        """

        self.report('SiestaCalculation<{}> did not reach geometry convergence'.format(node.pk))

        # We need to take care here of passing the
        # output geometry of old_calc to the new calculation
        if node.outputs.output_parameters.attributes["variable_geometry"]:
            self.ctx.inputs['structure'] = node.outputs.output_structure

        #The presence of `parent_calc_folder` triggers the real restart
        #meaning the copy of the .DM and the addition of `use-saved-dm` to the parameters
        self.ctx.inputs['parent_calc_folder'] = node.outputs.remote_folder

        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=80, exit_codes=_proc_exit_cod.SCF_NOT_CONV)  #pylint: disable = no-member
    def handle_error_scf_not_conv(self, node):
        """
        SCF convergence was not reached.  We need to restart from the
        previous calculation without changing any of the input parameters.
        """

        self.report('SiestaCalculation<{}> did not achieve scf convergence.'.format(node.pk))

        #We need to take care here of passing the
        #output geometry of old_calc to the new calculation
        if node.outputs.output_parameters.attributes["variable_geometry"]:
            self.ctx.inputs['structure'] = node.outputs.output_structure

        #The presence of `parent_calc_folder` triggers the real restart
        #meaning the copy of the .DM and the addition of use-saved-dm to the parameters
        self.ctx.inputs['parent_calc_folder'] = node.outputs.remote_folder

        #Should be also increase the number of scf max iterations?

        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=90, exit_codes=_proc_exit_cod.SPLIT_NORM)  #pylint: disable = no-member
    def handle_error_split_norm(self, node):
        """
        The split_norm parameter was too small.  We need to change it and restart.
        The minimum split_norm is stored in the logs of the old calculation.
        """

        from aiida_siesta.calculations.tkdict import FDFDict

        self.report('SiestaCalculation<{}> crashed with split_norm issue.'.format(node.pk))

        #This error happens only at the beginning of the run, therefore no real restart needed.
        #Just a new calculation with a new split_norm.
        #self.ctx.inputs['parent_calc_folder'] = node.outputs.remote_folder

        #Retrive the minimum split norm from the logs of failed calc.
        logs = orm.Log.objects.get_logs_for(node)
        for log in logs:
            if "Error in split_norm option" in log.message:
                mylog = log.message.split()
        new_split_norm = float(mylog[-1]) + 0.001

        #We want to understand the presence of "pao-split-norm" in input and:
        #1) if present, we change its value to the minimum allowed
        #2) if not present, we activate pao-SplitTailNorm
        #As we don't know in which sintax the user passed "pao-split-norm (remember
        #that every fdf variant is allowed), we translate the original dict in
        #a FDFDict that is aware of the equivalent keyword.
        #A FDFDict is not accepted in the context, but it is accepted in orm.Dict.
        transl_basis = FDFDict(self.ctx.inputs["basis"])
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
        self.ctx.inputs["basis"] = new_basis.get_dict()

        return ProcessHandlerReport(do_break=True)

    @process_handler(priority=89, exit_codes=_proc_exit_cod.BASIS_POLARIZ)  #pylint: disable = no-member
    def handle_error_basis_pol(self, node):  #pylint: disable = unused-argument
        """
        For the moment, we don't handle this error, but we terminate the WorkChain with
        a specific error code.
        """
        return ProcessHandlerReport(True, self.exit_codes.ERROR_BASIS_POL)

    #pylint: disable = no-member
    @process_handler(priority=60, exit_codes=[_proc_exit_cod.BANDS_PARSE_FAIL, _proc_exit_cod.BANDS_FILE_NOT_PRODUCED])
    def handle_error_bands(self, node):  #pylint: disable = unused-argument
        """
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
        from aiida_siesta.utils.inputs_generators import BaseWorkChainInputsGenerator
        return BaseWorkChainInputsGenerator(cls)
