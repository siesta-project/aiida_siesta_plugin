# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Bool, Int, Str, Float
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.remote import RemoteData

from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction
from aiida.common.links import LinkType

from aiida_siesta.data.psf import get_pseudos_from_structure

from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.calculations.stm import STMCalculation

                        
class SiestaSTMWorkChain(WorkChain):
    """
    STM Workchain. An example of workflow composition.
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
            cls.run_relax_and_analyze,
            cls.run_stm,   # We can run this directly, a combined scf+bands
            cls.run_results,
        )
        spec.dynamic_output()
                                         
    def setup_protocol(self):
        """
        Setup of context variables and inputs for the SiestaBaseWorkChain. Based on the specified
        protocol, we define values for variables that affect the execution of the calculations
        """
        self.ctx.inputs = {
            'code': self.inputs.code,
            'stm_code': self.inputs.stm_code,
            'height': self.inputs.height,
            'e1': self.inputs.e1,
            'e2': self.inputs.e2,
            'parameters': {},
            'settings': {},
            'options': ParameterData(dict={
                'resources': {
                    'num_machines': 1
                },
                'max_wallclock_seconds': 1800,
            }),
        }

        if self.inputs.protocol == 'standard':
            self.report('running the workchain in the "{}" protocol'.format(self.inputs.protocol.value))
            self.ctx.protocol = {
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.2,
                'dm_convergence_threshold': 1.0e-4,
                'forces_convergence_threshold': "0.02 eV/Ang",
                'min_meshcutoff': 100, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'md-type-of-run': "cg",
                'md-num-cg-steps': 10,
                'pseudo_familyname': 'lda-ag',
                # Future expansion. Add basis info, caveats, etc
                'atomic_heuristics': {
                    'H': { 'cutoff': 100 },
                    'Si': { 'cutoff': 100 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'DZP'
                }
                          
            }
        elif self.inputs.protocol == 'fast':
            self.report('running the workchain in the "{}" protocol'.format(self.inputs.protocol.value))
            self.ctx.protocol = {
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.25,
                'dm_convergence_threshold': 1.0e-3,
                'forces_convergence_threshold': "0.2 eV/Ang",
                'min_meshcutoff': 80, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'md-type-of-run': "cg",
                'md-num-cg-steps': 8,
                'pseudo_familyname': 'lda-ag',
                # Future expansion. Add basis info, caveats, etc
                'atomic_heuristics': {
                    'H': { 'cutoff': 50 },
                    'Si': { 'cutoff': 50 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'SZP'
                }
                          
            }
        else:
            self.abort_nowait('Protocol {} not known'.format(self.ctx.protocol.value))

    def setup_structure(self):
        """
        Very simple. Avoid seekpath for now
        """
        self.report('Running setup_structure')

        self.ctx.structure_initial_primitive = self.inputs.structure

    def setup_kpoints(self):
        """
        Define the k-point mesh for the relax calculation.
        """
        self.report('Running setup_kpoints')
        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(self.ctx.structure_initial_primitive)
        kpoints_mesh.set_kpoints_mesh_from_density(
            distance=self.ctx.protocol['kpoints_mesh_density'],
            offset=self.ctx.protocol['kpoints_mesh_offset']
        )
        
        self.ctx.kpoints_mesh = kpoints_mesh
        
    def setup_pseudo_potentials(self):
        """
        Based on the given input structure, get the 
        pseudo potentials for the different elements in the structure
        """
        self.report('Running setup_pseudo_potentials')
        structure = self.ctx.structure_initial_primitive
        pseudo_familyname = self.ctx.protocol['pseudo_familyname']
        self.ctx.inputs['pseudos'] = get_pseudos_from_structure(structure, pseudo_familyname)

    def setup_parameters(self):
        """
        Setup the default input parameters required for a
        SiestaCalculation and the SiestaBaseWorkChain
        """

        self.report('Running setup_parameters')
        structure = self.ctx.structure_initial_primitive
        meshcutoff = 0.0

        for kind in structure.get_kind_names():
            try:
                cutoff = self.ctx.protocol['atom_heuristics'][kind]['cutoff']
                meshcutoff = max(meshcutoff,cutoff)
            except:
                pass    # No problem. No heuristics, no info

        # In case we did not get anything, set a minimum value
        meshcutoff = max(self.ctx.protocol['min_meshcutoff'], meshcutoff)    
                
        self.ctx.inputs['parameters'] = {
            'dm-tolerance': self.ctx.protocol['dm_convergence_threshold'],
            'md-max-force-tol': self.ctx.protocol['forces_convergence_threshold'],
            'mesh-cutoff': "{} Ry".format(meshcutoff),
            'electronic-temperature': self.ctx.protocol['electronic_temperature'],
            'md-type-of-run': self.ctx.protocol['md-type-of-run'],
            'md-num-cg-steps': self.ctx.protocol['md-num-cg-steps']
        }

    def setup_basis(self):
        """
        Setup the basis dictionary.
        Very simple for now. Just the same for all elements. With more heuristics, we could do more.
        """
        self.report('Running setup_basis')
        self.ctx.inputs['basis'] = self.ctx.protocol['basis']
        
    def run_relax_and_analyze(self):
        """
        Run the SiestaBaseWorkChain to (relax) and analyze the input structure
        """
        self.report('Running run_relax_and_analyze')
        
        inputs = dict(self.ctx.inputs)
        
        ldos_e = "\n {e1} {e2} eV".format(e1=self.inputs.e1,e2=self.inputs.e2)
        inputs['parameters']['%block local-density-of-states'] = ldos_e


        # Final input preparation, wrapping dictionaries in ParameterData nodes
        # The code and options were set above
        # Pseudos was set above in ctx.inputs, and so in inputs
        
        inputs['kpoints'] = self.ctx.kpoints_mesh
        inputs['basis'] = ParameterData(dict=inputs['basis'])
        inputs['structure'] = self.ctx.structure_initial_primitive
        inputs['parameters'] = ParameterData(dict=inputs['parameters'])
        inputs['settings'] = ParameterData(dict=inputs['settings'])
        inputs['clean_workdir'] = Bool(False)
        inputs['max_iterations'] = Int(20)
        
        running = submit(SiestaBaseWorkChain, **inputs)
        self.report('launched SiestaBaseWorkChain<{}> in relax+ldos mode'.format(running.pid))
        
        return ToContext(workchain_relax=running)

    def run_stm(self):
        """
        Run a STMCalculation with the relaxed_calculation parent folder
        """
        self.report('Running stm calculation')
        # Get the remote folder of the last calculation in the previous workchain
        remote_folder = self.ctx.workchain_relax.get_outputs_dict()['remote_folder']

        stm_inputs = {}
        stm_inputs['code'] = self.ctx.inputs['stm_code']
        stm_inputs['parent_folder'] = remote_folder

        # Height of image plane, in Ang
        stm_inputs['parameters'] = ParameterData(dict={ 'z': self.ctx.inputs['height']})
        
        # Dummy dict
        settings_dict = {'a': 'b'}
        stm_inputs['settings'] = ParameterData(dict= settings_dict)
        
        # This should be just a dictionary!
        stm_inputs['_options'] = {
            'resources': {
                'num_machines': 1
            },
            'max_wallclock_seconds': 600
        }

        process = STMCalculation.process()
        running = submit(process, **stm_inputs)
        
        self.report('launching STMCalculation<{}>'.format(running.pid))
        
        return ToContext(stm_calc=running)

    def run_results(self):
        """
        Attach the relevant output nodes from the stm calculation to the workchain outputs
        for convenience
        """
        calculation_stm = self.ctx.stm_calc
        output_structure = self.ctx.workchain_relax.get_outputs_dict()['output_structure']

        self.report('workchain succesfully completed'.format())
        self.out('stm_array', calculation_stm.out.stm_array)
        self.out('output_structure', output_structure)
