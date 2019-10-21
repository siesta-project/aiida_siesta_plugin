# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida.orm import Code
from aiida.orm import (Int, Str, Float, Bool, Dict, StructureData, KpointsData,
                       ArrayData)
from aiida.engine import WorkChain, ToContext

from aiida_siesta.data.psf import get_pseudos_from_structure
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.calculations.stm import STMCalculation


class SiestaSTMWorkChain(WorkChain):
    """
    STM Workchain. An example of workflow composition.
    A separate workflow is only really needed if we desire to have
    a robust scf and/or relaxation cycle.
    """
    def __init__(self, *args, **kwargs):
        super(SiestaSTMWorkChain, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(SiestaSTMWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('stm_code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input('protocol', valid_type=Str, default=Str('standard'))
        spec.input('height', valid_type=Float)
        spec.input('e1', valid_type=Float)
        spec.input('e2', valid_type=Float)
        spec.outline(
            cls.setup_protocol,
            cls.setup_structure,
            cls.setup_kpoints,
            cls.setup_pseudo_potentials,
            cls.setup_parameters,
            cls.setup_basis,
            cls.run_relax,  # This should be optional
            cls.
            run_scf_and_ldos,  # We could run this directly, a combined scf+ldos
            cls.run_stm,
            cls.run_results,
        )

        spec.output('stm_array', valid_type=ArrayData)
        # These will be inherited from the Base workchain output
        #spec.output('output_structure',
        #            valid_type=StructureData,
        #            required=False)

        # This could be 'output_parameters', also inherited from the Base workchain, representing
        # the summary of energies, etc. coming out of the last SiestaCalculation run.
        # but its relevance is limited for a "stm" workflow. Left here for now for debugging.
        #spec.output('output_parameters', valid_type=Dict)

        spec.exit_code(140,
                       'ERROR_PROTOCOL_NOT_FOUND',
                       message='The protocol specified is not known')
        spec.exit_code(
            160,
            'ERROR_RELAXED_STRUCTURE_NOT_AVAILABLE',
            message='Failed to get the output structure from the relaxation run'
        )

    def setup_protocol(self):
        """
        Setup of context variables and inputs for the SistaBaseWorkChain. Based on the specified
        protocol, we define values for variables that affect the execution of the calculations
        """
        self.ctx.inputs = {
            'code': self.inputs.code,
            'parameters': {},
            'settings': {},
            'options': {
                'resources': {
                    'num_machines': 1
                },
                'max_wallclock_seconds': 1800,
            },
        }
        #
        #  Separate inputs for LDOS generation and STM plot generation
        #
        self.ctx.inputs_stm = {
            'stm_code': self.inputs.stm_code,
            'height': self.inputs.height,
            'e1': self.inputs.e1,
            'e2': self.inputs.e2,
            'options': {
                'resources': {
                    'num_machines': 1
                },
                'max_wallclock_seconds': 1800,
            }
        }

        if self.inputs.protocol == 'standard':
            self.report('Using the "{}" protocol'.format(
                self.inputs.protocol.value))
            self.ctx.protocol = {
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.2,
                'dm_convergence_threshold': 1.0e-4,
                'forces_convergence_threshold': "0.02 eV/Ang",
                'min_meshcutoff': 100,  # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'md-type-of-run': "cg",
                'md-num-cg-steps': 10,
                'pseudo_familyname': 'sample_psf_family',
                # Future expansion. Add basis info, caveats, etc
                'atomic_heuristics': {
                    'H': {
                        'cutoff': 100
                    },
                    'Si': {
                        'cutoff': 100
                    }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'DZP'
                }
            }

        elif self.inputs.protocol == 'fast':
            self.report('Using the "{}" protocol'.format(
                self.inputs.protocol.value))
            self.ctx.protocol = {
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.25,
                'dm_convergence_threshold': 1.0e-3,
                'forces_convergence_threshold': "0.2 eV/Ang",
                'min_meshcutoff': 80,  # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'md-type-of-run': "cg",
                'md-num-cg-steps': 8,
                'pseudo_familyname': 'sample_psf_family',
                # Future expansion. Add basis info, caveats, etc
                'atomic_heuristics': {
                    'H': {
                        'cutoff': 50
                    },
                    'Si': {
                        'cutoff': 50
                    }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'SZP'
                }
            }
        else:
            self.report('Protocol {} not known'.format(
                self.inputs.protocol.value))
            return self.exit_codes.ERROR_PROTOCOL_NOT_FOUND

    def setup_structure(self):
        """
        Just a stub for future expansion, maybe with normalization
        """

        self.ctx.structure_initial_primitive = self.inputs.structure

    def setup_kpoints(self):
        """
        Define the k-point mesh for the relax and scf calculations.
        """

        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(
            self.ctx.structure_initial_primitive)
        kpoints_mesh.set_kpoints_mesh_from_density(
            distance=self.ctx.protocol['kpoints_mesh_density'],
            offset=self.ctx.protocol['kpoints_mesh_offset'])

        self.ctx.kpoints_mesh = kpoints_mesh

    def setup_pseudo_potentials(self):
        """
        Based on the given input structure and the protocol, use the SSSP library to determine the
        optimal pseudo potentials for the different elements in the structure
        """
        structure = self.ctx.structure_initial_primitive
        pseudo_familyname = self.ctx.protocol['pseudo_familyname']
        self.ctx.inputs['pseudos'] = get_pseudos_from_structure(
            structure, pseudo_familyname)

    def setup_parameters(self):
        """
        Setup the default input parameters required for a SiestaCalculation and the SiestaBaseWorkChain
        """
        structure = self.ctx.structure_initial_primitive
        meshcutoff = 0.0

        for kind in structure.get_kind_names():
            try:
                cutoff = self.ctx.protocol['atom_heuristics'][kind]['cutoff']
                meshcutoff = max(meshcutoff, cutoff)
            except:
                pass  # No problem. No heuristics, no info

        meshcutoff = max(
            self.ctx.protocol['min_meshcutoff'],
            meshcutoff)  # In case we did not get anything, set a minimum value

        self.ctx.inputs['parameters'] = {
            'dm-tolerance': self.ctx.protocol['dm_convergence_threshold'],
            'md-max-force-tol':
            self.ctx.protocol['forces_convergence_threshold'],
            'mesh-cutoff': "{} Ry".format(meshcutoff),
            'electronic-temperature':
            self.ctx.protocol['electronic_temperature'],
            'md-type-of-run': self.ctx.protocol['md-type-of-run'],
            'md-num-cg-steps': self.ctx.protocol['md-num-cg-steps']
        }

    def setup_basis(self):
        """
        Setup the basis dictionary.
        Very simple for now. Just the same for all elements. With more heuristics, we could do more.
        """
        self.ctx.inputs['basis'] = self.ctx.protocol['basis']

    def run_relax(self):
        """
        Run the SiestaBaseWorkChain to relax the input structure
        """

        inputs = dict(self.ctx.inputs)

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        # The code and options (_options?)  were set above
        # Pseudos was set above in 'ctx.inputs', and so it is in 'inputs' already

        inputs['kpoints'] = self.ctx.kpoints_mesh
        inputs['basis'] = Dict(dict=inputs['basis'])
        inputs['structure'] = self.ctx.structure_initial_primitive
        inputs['parameters'] = Dict(dict=inputs['parameters'])
        inputs['settings'] = Dict(dict=inputs['settings'])
        inputs['options'] = Dict(dict=inputs['options'])
        inputs['clean_workdir'] = Bool(False)
        inputs['max_iterations'] = Int(20)

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report(
            'launched SiestaBaseWorkChain<{}> in relaxation mode'.format(
                running.pk))

        return ToContext(workchain_relax=running)

    def run_scf_and_ldos(self):
        """
        Run the SiestaBaseWorkChain in scf+ldos mode on the primitive cell of the relaxed input structure
        """

        try:
            structure = self.ctx.workchain_relax.outputs.output_structure
        except:
            return self.exit_codes.ERROR_RELAXED_STRUCTURE_NOT_AVAILABLE

        # Do we need further refinement by Seekpath on this=? (eventually)?
        self.ctx.structure_relaxed_primitive = structure

        inputs = dict(self.ctx.inputs)
        ldos_e = "\n {e1} {e2} eV \n %endblock local-density-of-states".format(
            e1=self.inputs.e1.value, e2=self.inputs.e2.value)
        inputs['parameters']['%block local-density-of-states'] = ldos_e

        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(
            self.ctx.structure_relaxed_primitive)
        kpoints_mesh.set_kpoints_mesh_from_density(
            distance=self.ctx.protocol['kpoints_mesh_density'],
            offset=self.ctx.protocol['kpoints_mesh_offset'])

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        inputs['kpoints'] = kpoints_mesh
        inputs['structure'] = self.ctx.structure_relaxed_primitive
        inputs['parameters'] = Dict(dict=inputs['parameters'])
        inputs['basis'] = Dict(dict=inputs['basis'])
        inputs['settings'] = Dict(dict=inputs['settings'])
        inputs['options'] = Dict(dict=inputs['options'])

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report('launched SiestaBaseWorkChain<{}> in scf+ldos mode'.format(
            running.pk))

        return ToContext(workchain_base_ldos=running)

    def run_stm(self):
        """
        Run a STMCalculation with the relaxed_calculation parent folder
        """
        self.report('Running stm calculation')

        # Get the remote folder of the last calculation in the previous workchain
        base_ldos = self.ctx.workchain_base_ldos
        remote_folder = base_ldos.outputs.remote_folder

        stm_inputs = {}
        stm_inputs['code'] = self.ctx.inputs_stm['stm_code']
        stm_inputs['ldos_folder'] = remote_folder

        # Height of image plane, in Ang
        stm_inputs['parameters'] = Dict(
            dict={'z': self.ctx.inputs_stm['height']})

        # The "metadata.options.resources" seems to be needed
        metadata = {"options": self.ctx.inputs_stm['options']}
        stm_inputs['metadata'] = metadata  # No Dict... ??

        # Dummy dict
        settings_dict = {'a': 'b'}
        stm_inputs['settings'] = Dict(dict=settings_dict)

        running = self.submit(STMCalculation, **stm_inputs)
        self.report('launching STMCalculation<{}>'.format(running.pk))

        return ToContext(stm_calc=running)

    def run_results(self):
        """
        Attach the relevant output nodes
        """

        from aiida.engine import ExitCode

        stm_plot_calc = self.ctx.stm_calc
        stm_array = stm_plot_calc.outputs.stm_array
        self.out('stm_array', stm_array)

        # If we end up keeping this workflow as a relax+stm_image task, we should get the relaxed
        # structure.
        # output_structure = self.ctx.workchain_relax.get_outputs_dict()['output_structure']
        #self.out('output_structure', output_structure)

        self.report('STM workchain succesfully completed')
        return ExitCode(0)
