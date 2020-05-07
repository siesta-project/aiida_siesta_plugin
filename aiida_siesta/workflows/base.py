from aiida import orm
from aiida.engine import BaseRestartWorkChain, ProcessHandlerReport, process_handler, while_
from aiida_siesta.data.common import get_pseudos_from_structure
from aiida_siesta.calculations.siesta import SiestaCalculation


def prepare_pseudo_inputs(structure, pseudos, pseudo_family):

    if pseudos and pseudo_family:
        raise ValueError('you cannot specify both "pseudos" and "pseudo_family"')
    elif pseudos is None and pseudo_family is None:
        raise ValueError('neither an explicit pseudos dictionary nor a pseudo_family was specified')
    elif pseudo_family:
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
        #Required by the plugin
        spec.input('parameters', valid_type=orm.Dict)
        spec.input('basis', valid_type=orm.Dict, required=False)
        spec.input('settings', valid_type=orm.Dict, required=False)
        #Required by any CalcJob
        spec.input('options', valid_type=orm.Dict)

        spec.outline(
            cls.setup,
            cls.prepare_inputs,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )

        spec.output('forces_and_stress', valid_type=orm.ArrayData, required=False)
        spec.output('bands', valid_type=orm.BandsData, required=False)
        spec.output('output_structure', valid_type=orm.StructureData, required=False)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('remote_folder', valid_type=orm.RemoteData)

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

    @process_handler(exit_codes=SiestaCalculation.exit_codes.GEOM_NOT_CONV)  #pylint: disable = no-member
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

    @process_handler(exit_codes=SiestaCalculation.exit_codes.SCF_NOT_CONV)  #pylint: disable = no-member
    def handle_error_scf_not_conv(self, node):
        """
        SCF convergence was not reached.  We need to restart from the
        previous calculation without changing any of the input parameters.
        """

        self.report('SiestaCalculation<{}> did not achieve scf convergence.'.format(node.pk))

        #The presence of `parent_calc_folder` triggers the real restart
        #meaning the copy of the .DM and the addition of use-saved-dm to the parameters
        self.ctx.inputs['parent_calc_folder'] = node.outputs.remote_folder

        return ProcessHandlerReport(do_break=True)
