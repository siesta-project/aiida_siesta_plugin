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
from aiida_siesta.calculations.fcbuild import FCBuildCalculation
from aiida_siesta.calculations.vibrator import VibratorCalculation

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.12.0"
__contributors__ = "Victor M. Garcia-Suarez, ..."
                        
class SiestaVibraWorkChain(WorkChain):
    """
    Vibra Workchain. An example of workflow composition for phonons
    """

    def __init__(self, *args, **kwargs):
        super(SiestaVibraWorkChain, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(SiestaVibraWorkChain, cls).define(spec)
        spec.input('fcbuild_code', valid_type=Code)
        spec.input('code', valid_type=Code)
        spec.input('vibrator_code', valid_type=Code)
        spec.input('fcbparams', valid_type=ParameterData)
        spec.input('structure', valid_type=StructureData)
        spec.input('protocol', valid_type=Str, default=Str('standard'))
        spec.input('bandskpoints', valid_type=KpointsData)
        spec.outline(
            cls.setup_structure,
            cls.run_fcbuild,
            cls.setup_rra_inputs,
            cls.setup_protocol,
            cls.setup_pseudo_potentials,
            cls.setup_parameters,
            cls.setup_basis,
            cls.setup_kpoints,
            cls.run_relax_and_analyze,
            cls.run_vibrator,
            cls.run_results,
        )
        spec.dynamic_output()
                                         
    def setup_structure(self):
        """
        Very simple. Avoid seekpath for now
        """
        self.report('Running setup_structure')

        self.ctx.structure_initial_primitive = self.inputs.structure

    def run_fcbuild(self):
        """
        Run an initial FCBuildCalculation
        """
        self.report('Running fcbuild calculation')

        fcbuild_inputs = {}
        fcbuild_inputs['code'] = self.inputs.fcbuild_code
        fcbuild_inputs['structure'] = self.ctx.structure_initial_primitive
        fcbuild_inputs['parameters'] = self.inputs.fcbparams

        settings_dict = {}
        fcbuild_inputs['settings'] = ParameterData(dict=settings_dict)

        fcbuild_inputs['_options'] = {
            'resources': {
                'num_machines': 1
            },
            'max_wallclock_seconds': 600
        }
        
        process = FCBuildCalculation.process()
        running = submit(process, **fcbuild_inputs)
        
        self.report('launching FCBuildCalculation<{}>'.format(running.pid))
        
        return ToContext(fcbuild_calc=running)

    def setup_rra_inputs(self):
        """
        Setup of context variables and inputs for the SiestaBaseWorkChain. Based on the specified
        protocol, we define values for variables that affect the execution of the calculations
        """
        self.ctx.rra_inputs = {
            'code': self.inputs.code,
            'parameters': {},
            'settings': {},
            'options': ParameterData(dict={
                'resources': {
                    'num_machines': 1
                },
                'max_wallclock_seconds': 1800,
            }),
        }

    def setup_protocol(self):
        if self.inputs.protocol == 'standard':
            self.report('running the workchain in the "{}" protocol'.format(self.inputs.protocol.value))
            self.ctx.protocol = {
                'kpoints_mesh': 1,
                'dm_convergence_threshold': 1.0e-4,
                'min_meshcutoff': 100, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'pseudo_familyname': 'lda-ag',
                'atomic_heuristics': {
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
                'kpoints_mesh': 1,
                'dm_convergence_threshold': 1.0e-3,
                'min_meshcutoff': 80, # In Rydberg (!)
                'electronic_temperature': "25.0 meV",
                'pseudo_familyname': 'lda-ag',
                'atomic_heuristics': {
                    'Si': { 'cutoff': 50 }
                },
                'basis': {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'SZP'
                }
                          
            }
        else:
            self.abort_nowait('Protocol {} not known'.format(self.ctx.protocol.value))

    def setup_pseudo_potentials(self):
        """
        Based on the given input structure, get the 
        pseudo potentials for the different elements in the structure
        """
        self.report('Running setup_pseudo_potentials')
        structure = self.ctx.structure_initial_primitive
        pseudo_familyname = self.ctx.protocol['pseudo_familyname']
        self.ctx.rra_inputs['pseudos'] = get_pseudos_from_structure(structure, pseudo_familyname)

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
                
        self.ctx.rra_inputs['parameters'] = {
            'dm-tolerance': self.ctx.protocol['dm_convergence_threshold'],
            'mesh-cutoff': "{} Ry".format(meshcutoff),
            'electronic-temperature': self.ctx.protocol['electronic_temperature'],
            '%include': "FC.fdf"
        }

    def setup_basis(self):
        """
        Setup the basis dictionary.
        Very simple for now. Just the same for all elements. With more heuristics, we could do more.
        """
        self.report('Running setup_basis')
        self.ctx.rra_inputs['basis'] = self.ctx.protocol['basis']
        
    def setup_kpoints(self):
        """
        Define the k-point mesh for the relax calculation.
        """
        self.report('Running setup_kpoints')
        kpoints_mesh = KpointsData()
        kpmesh=self.ctx.protocol['kpoints_mesh']
        kpoints_mesh.set_kpoints_mesh([kpmesh,kpmesh,kpmesh])
        
        self.ctx.kpoints_mesh = kpoints_mesh

    def run_relax_and_analyze(self):
        """
        Run the SiestaBaseWorkChain to (relax) and analyze the input structure
        """
        self.report('Running run_relax_and_analyze')

        rra_inputs = {}
        rra_inputs = dict(self.ctx.rra_inputs)

        # Get the remote folder of the last calculation in the previous workchain
        remote_folder = self.ctx.fcbuild_calc.out.remote_folder
        rra_inputs['parent_folder'] = remote_folder

        rra_inputs['kpoints'] = self.ctx.kpoints_mesh
        rra_inputs['basis'] = ParameterData(dict=rra_inputs['basis'])
        rra_inputs['structure'] = self.ctx.structure_initial_primitive
        rra_inputs['parameters'] = ParameterData(dict=rra_inputs['parameters'])
        rra_inputs['settings'] = ParameterData(dict=rra_inputs['settings'])
        #rra_inputs['settings']['ADDITIONAL_REMOTE_COPY_LIST'] = [('aiida.FC', './')]
        rra_inputs['clean_workdir'] = Bool(False)
        rra_inputs['max_iterations'] = Int(20)
        
        running = submit(SiestaBaseWorkChain, **rra_inputs)
        self.report('launched SiestaBaseWorkChain<{}> in relax+analyze mode'.format(running.pid))
        
        return ToContext(workchain_relax=running)

    def run_vibrator(self):
        """
        Run a VibratorCalculation with the calculation parent folder
        """
        self.report('Running vibrator calculation')
        # Get the remote folder of the last calculation in the previous workchain
        #remote_folder = self.ctx.workchain_relax.get_outputs_dict()['remote_folder']
        remote_folder = self.ctx.workchain_relax.out.remote_folder
        #remote_folder = self.ctx.fcbuild_calc.out.remote_folder

        vibrator_inputs = {}
        vibrator_inputs['code'] = self.inputs.vibrator_code
        vibrator_inputs['parent_folder'] = remote_folder
        vibrator_inputs['structure'] = self.inputs.structure
        vibrator_inputs['bandskpoints'] = self.inputs.bandskpoints

        settings_dict = {}
        vibrator_inputs['settings'] = ParameterData(dict=settings_dict)

        parameters_dict = {}
        vibrator_inputs['parameters'] = ParameterData(dict=parameters_dict)

        vibrator_inputs['_options'] = {
            'resources': {
                'parallel_env': 'mpi',
                'tot_num_mpiprocs':1
            },
            'max_wallclock_seconds': 600
        }

        process = VibratorCalculation.process()
        running = submit(process, **vibrator_inputs)
        
        self.report('launching VibratorCalculation<{}>'.format(running.pid))
        
        return ToContext(vibrator_calc=running)

    def run_results(self):
        """
        Attach the relevant output nodes from the vibra calculation to the workchain outputs
        for convenience
        """
        calculation_vibra = self.ctx.vibrator_calc
        #calculation_vibra = self.ctx.fcbuild_calc
        #output_structure = self.ctx.fcbuild_calc.get_outputs_dict()['output_structure']

        self.report('workchain succesfully completed'.format())
        #self.out('fcbuild_array', calculation_vibra.out.fcbuild_array)
        #self.out('output_structure', output_structure)
