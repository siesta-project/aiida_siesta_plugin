from aiida import orm
from aiida.engine import WorkChain, calcfunction, ToContext
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.orm.nodes.data.structure import Kind

@calcfunction
def parse_neb(retrieved,ref_structure):
    """
    Wrapper to preserve provenance.
    :param: retrieved:  the retrieved folder from a NEB calculation
                        (containing .xyz files and NEB data files)
    :param: ref_structure: a reference structure
    :return: a Trajectory object generated from the .xyz files, and
             with extra arrays for NEB results.
    """
    import os
    from aiida.orm import TrajectoryData
    from aiida_siesta.utils.xyz_utils import get_structure_list_from_folder
    from aiida_siesta.utils.neb import parse_neb_results
    
    folder_path = retrieved._repository._get_base_folder().abspath
    struct_list = get_structure_list_from_folder(folder_path,ref_structure)

    traj = TrajectoryData(struct_list)

    neb_results_file = 'NEB.results'
    if neb_results_file in retrieved._repository.list_object_names():
        neb_results_path = os.path.join(folder_path, neb_results_file)
        annotated_traj = parse_neb_results(neb_results_path,traj)

        _kinds_raw = [ k.get_raw() for k in ref_structure.kinds ]
        annotated_traj.set_attribute('kinds', _kinds_raw)

    return annotated_traj

class SiestaBaseNEBWorkChain(WorkChain):

    """
    Workchain to run a NEB MEP optimization
    starting from a guessed path
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(SiestaCalculation,
                           exclude=('structure','lua','metadata'))

        # We might enforce the kinds annotation by using a new data type,
        # but it would be wasteful
        spec.input('starting_path', valid_type=orm.TrajectoryData,
                   help='Starting Path')
        spec.input('neb_script', valid_type=orm.SinglefileData,
                   help='Lua script for NEB engine')

        spec.input('options', valid_type=orm.Dict, help='Options')

        # These, together with n_images, are passed as 'lua' parameters
        
        spec.input('spring_constant', valid_type=orm.Float, required=False)
        # spec.input('climbing_image', valid_type=orm.Bool, required=False)
        # spec.input('max_number_of_neb_iterations', valid_type=orm.Int, required=False)
        # ... tolerances, etc are encoded in the Siesta params dictionary.
        
        
        spec.output('neb_output_package', valid_type=orm.TrajectoryData)

        spec.outline(

            cls.check_input_path,
            cls.run_neb,
            cls.run_results,
        )
        spec.exit_code(201, 'ERROR_PATH_SPEC', message='The path specification is faulty')
        spec.exit_code(201, 'ERROR_NEB_CALC', message='The NEB calculation failed')

    def check_input_path(self):
        """
        Make sure that the input set of images is consistent, and is annotated with
        the kinds of the structure
        """
        path = self.inputs.starting_path

        if path.numsteps == 0:
            self.report('The trajectory data object does not contain any structures')
            return self.exit_codes.ERROR_PATH_SPEC

        if path.numsteps == 1:
            self.report('The trajectory data object does not represent a path')
            return self.exit_codes.ERROR_PATH_SPEC
            
        if path.numsteps == 2:
            self.report('The trajectory data object contains only two structures...')
            # We could just interpolate, but here this is an error
            return self.exit_codes.ERROR_PATH_SPEC

        # ... further "smoothness" tests could be implemented if needed
        
        try:
            _kinds_raw = path.get_attribute('kinds')
        except AttributeError:
            self.report('No kinds attribute found in TrajectoryData object')
            return self.exit_codes.ERROR_PATH_SPEC

        # Create proper kinds list from list of raw dictionaries
        _kinds = [ Kind(raw=kr) for kr in _kinds_raw]
        
        ref_structure = path.get_step_structure(0,custom_kinds=_kinds)
        
        self.ctx.reference_structure = ref_structure
        
    def run_neb(self):
        """
        Run a SiestaCalculation with a specific NEB images
        input.
        """
        
        inputs = self.exposed_inputs(SiestaCalculation)

        neb_path = self.inputs.starting_path
        kinds = self.ctx.reference_structure.kinds
        ghost_dict = self.inputs.basis
        neb_image_prefix = 'image-'
            
        from aiida.orm import FolderData
        from aiida.common.folders import SandboxFolder
        from aiida_siesta.utils.xyz_utils import write_xyz_file_from_structure
        # Temporary folder
        with SandboxFolder() as folder:
            folder_path= folder.abspath

            # loop over structures and create xyz files
            for i in range(neb_path.numsteps):
                s_image_tmp = neb_path.get_step_structure(i,custom_kinds=kinds)
                s_image, dummy = add_ghost_sites_to_structure(s_image_tmp, ghost_dict)
                filename= folder.get_abs_path("{}{}.xyz".format(neb_image_prefix,i))
                write_xyz_file_from_structure(s_image,filename)

            lua_input_files=FolderData(tree=folder_path)

        n_images = self.inputs.starting_path.numsteps - 2
        spring_constant = 0.1
        if 'spring_constant' in self.inputs:
            spring_constant = self.inputs.spring_constant.value

        lua_parameters = { }
        lua_parameters['neb_spring_constant'] = spring_constant
        lua_parameters['number_of_internal_images_in_path'] = n_images
        lua_parameters['neb_image_file_prefix'] = 'image-'
        
        inputs['lua'] = {}
        inputs['lua']['script'] = self.inputs.neb_script
        inputs['lua']['input_files'] = lua_input_files
        inputs['lua']['parameters'] = orm.Dict(dict=lua_parameters)
        inputs['lua']['retrieve_list'] = orm.List(list=['*.xyz','NEB.results'])
        
        inputs['structure'] = self.ctx.reference_structure

        #
        # Note this
        #
        inputs['metadata'] = {
            "label": "NEB calculation",
            'options': self.inputs.options.get_dict(),
         }

        running = self.submit(SiestaCalculation, **inputs)
        self.report(f'Launched SiestaCalculation<{running.pk}> for NEB.')

        return ToContext(neb_wk=running)

    def run_results(self):

        if not self.ctx.neb_wk.is_finished_ok:
                return self.exit_codes.ERROR_NEB_CALC
        else:
            outps = self.ctx.neb_wk.outputs

        # Use the 'retrieved' folder and parse the NEB data here
        #
        retrieved_folder = outps['retrieved']
        ref_structure = self.ctx.reference_structure

        annotated_traj = parse_neb(retrieved_folder, ref_structure)
        # Get better diagnostics from calcfunction...
        if annotated_traj.numsteps == 0:
            self.logger.warning("No .xyz files retrieved")
            return self.exit_codes.NO_NEB_XYZ_FILES

        self.out('neb_output_package', annotated_traj)
        n_iterations = annotated_traj.get_attribute('neb_iterations')
        self.report(f'NEB process done in {n_iterations} iterations.')
                

#    @classmethod
#    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
#        from aiida_siesta.utils.inputs_generators import BaseWorkChainInputsGenerator
#        return BaseWorkChainInputsGenerator(cls)
