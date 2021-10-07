from aiida import orm
from aiida.engine import WorkChain, ToContext
from aiida_siesta.workflows.neb_base import SiestaBaseNEBWorkChain
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.utils.structures import exchange_sites_in_structure
from aiida_siesta.utils.structures import compute_mid_path_position
from aiida_siesta.utils.structures import find_intermediate_structure
from aiida_siesta.utils.interpol import interpolate_two_structures_ase


class ExchangeBarrierWorkChain(WorkChain):
    """
    Workchain to compute the barrier for exchange of two atoms
    in a structure.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(SiestaBaseWorkChain, exclude=('structure',), namespace="initial")
        spec.expose_inputs(SiestaBaseWorkChain, exclude=('structure',), namespace="final")

        spec.expose_inputs(SiestaBaseNEBWorkChain, exclude=('starting_path',), namespace="neb")

        spec.input('initial_structure', valid_type=orm.StructureData, help='Initial structure')

        spec.input('first_index', valid_type=orm.Int, help='Index of first atom in structure')
        spec.input('second_index', valid_type=orm.Int, help='Index of second atom structure')
        spec.input('migration_direction', valid_type=orm.List, help='Migration direction (in lattice coordinates)')

        spec.input('n_images', valid_type=orm.Int, help='Number of (internal) images in Path (odd!!)')  # validate

        spec.expose_outputs(SiestaBaseNEBWorkChain)

        spec.outline(
            cls.prepare_structures, cls.relax_initial, cls.relax_final, cls.prepare_initial_path, cls.run_NEB_workchain,
            cls.check_results
        )

        spec.exit_code(200, 'ERROR_MAIN_WC', message='The end-point relaxation SiestaBaseWorkChain failed')
        spec.exit_code(250, 'ERROR_CONFIG', message='Cannot generate initial path correctly')
        spec.exit_code(300, 'ERROR_NEB_WK', message='SiestaBaseNEBWorkChain did not finish correctly')

    def prepare_structures(self):
        """
        Generate exchanged structure as final end-point
        """

        s_initial = self.inputs.initial_structure

        i1 = self.inputs.first_index.value
        i2 = self.inputs.second_index.value

        s_final = exchange_sites_in_structure(s_initial, i1, i2)

        self.ctx.s_initial = s_initial
        self.ctx.s_final = s_final

        self.report(f'Created initial and final structures')

    def relax_initial(self):
        """
        Run the SiestaBaseWorkChain, might be a relaxation or a scf only.
        """

        inputs = self.exposed_inputs(SiestaBaseWorkChain, namespace='initial')
        inputs['structure'] = self.ctx.s_initial

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report(f'Launched SiestaBaseWorkChain<{running.pk}> to relax the initial structure.')

        return ToContext(initial_relaxation_wk=running)

    def relax_final(self):
        """
        Run the SiestaBaseWorkChain, might be a relaxation or a scf only.
        """

        inputs = self.exposed_inputs(SiestaBaseWorkChain, namespace='final')
        inputs['structure'] = self.ctx.s_final

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report(f'Launched SiestaBaseWorkChain<{running.pk}> to relax the final structure.')

        return ToContext(final_relaxation_wk=running)

    def prepare_initial_path(self):
        """
        Perhaps more heuristics are needed?
        Here we just interpolate.
        """

        initial_wk = self.ctx.initial_relaxation_wk
        if not initial_wk.is_finished_ok:
            return self.exit_codes.ERROR_MAIN_WC

        final_wk = self.ctx.final_relaxation_wk
        if not final_wk.is_finished_ok:
            return self.exit_codes.ERROR_MAIN_WC

        s_initial = initial_wk.outputs.output_structure
        s_final = final_wk.outputs.output_structure

        n_images = self.inputs.n_images.value

        #
        # Add here any heuristics, before handling the
        # path for further refinement
        #
        #---------------------------------------------
        # The basic heuristic here is to avoid head-on collissions
        # by defining an "avoidance cylinder" around the line
        # joining the two atoms exchanged. The input "migration_direction"
        # serves to define a point on the surface of that cylinder, at
        # the mid-point, which is used as the mid-point of the trial path.

        migration_direction = self.inputs.migration_direction.get_list()
        i1 = self.inputs.first_index.value
        i2 = self.inputs.second_index.value

        i1_mid_path_position = compute_mid_path_position(s_initial, i1, i2, migration_direction)

        #
        s_intermediate = find_intermediate_structure(s_initial, i1, i1_mid_path_position, i2)

        #
        # The starting_path is now built from two sections
        # We assume that the number of internal images is odd,
        # so that n_images // 2 is the number of internal images
        # of each section

        first_list = interpolate_two_structures_ase(s_initial, s_intermediate, n_images // 2)
        second_list = interpolate_two_structures_ase(s_intermediate, s_final, n_images // 2)

        #
        # Remove duplicate central point
        #
        images_list = first_list[:-1] + second_list

        if len(images_list) != n_images + 2:
            self.report(f"Number of images: {n_images} /= list length")
            return self.exit_codes.ERROR_CONFIG

        #
        # We might need a more general refiner, starting
        # with the trial path
        #
        # refined_path = refine_neb_path(starting_path)

        path_object = orm.TrajectoryData(images_list)
        #
        # Use a 'serializable' dictionary instead of the
        # actual kinds list
        #
        _kinds_raw = [k.get_raw() for k in s_initial.kinds]
        path_object.set_attribute('kinds', _kinds_raw)

        self.ctx.path = path_object

        self.report(f'Generated starting path for NEB.')

    def run_NEB_workchain(self):

        inputs = self.exposed_inputs(SiestaBaseNEBWorkChain, namespace='neb')

        print(inputs)

        inputs['starting_path'] = self.ctx.path

        running = self.submit(SiestaBaseNEBWorkChain, **inputs)

        self.report(f'Launched SiestaBaseNEBWorkChain<{running.pk}> to find MEP for atom exchange.')

        return ToContext(neb_wk=running)

    def check_results(self):
        """
        All checks are done in the NEB workchain
        """

        if not self.ctx.neb_wk.is_finished_ok:
            return self.exit_codes.ERROR_NEB_WK

        outps = self.ctx.neb_wk.outputs
        self.out('neb_output_package', outps['neb_output_package'])

        self.report(f'ExchangeBarrier workchain done.')


#    @classmethod
#    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
#        from aiida_siesta.utils.inputs_generators import BaseWorkChainInputsGenerator
#        return BaseWorkChainInputsGenerator(cls)
