# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Bool, Int, Str
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.work.run import submit, run
from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction
from aiida.common.links import LinkType

from aiida_siesta.data.psf import get_pseudos_from_structure
##from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida.orm.data.array.kpoints import KpointsData

class SiestaBandsWorkChain(WorkChain):
    """
    Bands Workchain. An example of workflow composition.
    A separate bands workflow is only needed if we desire to separate the
    relaxation from the final run.
    """

    def __init__(self, *args, **kwargs):
        super(SiestaBandsWorkChain, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(SiestaBandsWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input('protocol', valid_type=Str, default=Str('standard'))
        spec.outline(
            cls.setup_protocol,
            cls.setup_structure,
            cls.setup_kpoints,
            cls.setup_pseudo_potentials,
            cls.setup_parameters,
            cls.setup_basis,
            cls.run_relax,
            cls.run_bands,   # We can run this directly, a combined scf+bands
            cls.run_results,
        )
        spec.dynamic_output()
                                         
    def setup_protocol(self):
        """
        Setup of context variables and inputs for the SistaBaseWorkChain. Based on the specified
        protocol, we define values for variables that affect the execution of the calculations
        """
        self.ctx.inputs = {
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
        Just a stub for future expansion, maybe with normalization
        """

        self.ctx.structure_initial_primitive = self.inputs.structure

    def setup_kpoints(self):
        """
        Define the k-point mesh for the relax and scf calculations.
        """

        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(self.ctx.structure_initial_primitive)
        kpoints_mesh.set_kpoints_mesh_from_density(
            distance=self.ctx.protocol['kpoints_mesh_density'],
            offset=self.ctx.protocol['kpoints_mesh_offset']
        )
        
        self.ctx.kpoints_mesh = kpoints_mesh
        
    def setup_pseudo_potentials(self):
        """
        Based on the given input structure and the protocol, use the SSSP library to determine the
        optimal pseudo potentials for the different elements in the structure
        """
        structure = self.ctx.structure_initial_primitive
        pseudo_familyname = self.ctx.protocol['pseudo_familyname']
        self.ctx.inputs['pseudos'] = get_pseudos_from_structure(structure, pseudo_familyname)

    def setup_parameters(self):
        """
        Setup the default input parameters required for a SiestaCalculation and the SiestaBaseWorkChain
        """
        structure = self.ctx.structure_initial_primitive
        meshcutoff = 0.0

        for kind in structure.get_kind_names():
            try:
                cutoff = self.ctx.protocol['atom_heuristics'][kind]['cutoff']
                meshcutoff = max(meshcutoff,cutoff)
            except:
                pass    # No problem. No heuristics, no info

        meshcutoff = max(self.ctx.protocol['min_meshcutoff'], meshcutoff)    # In case we did not get anything, set a minimum value
                
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
        self.ctx.inputs['basis'] = self.ctx.protocol['basis']
        
    def run_relax(self):
        """
        Run the SiestaBaseWorkChain to relax the input structure
        """
        self.report('Running run_relax')
        inputs = dict(self.ctx.inputs)

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        # The code and options (_options?)  were set above
        # Pseudos was set above in 'ctx.inputs', and so it is in 'inputs' already
        
        inputs['kpoints'] = self.ctx.kpoints_mesh
        inputs['basis'] = ParameterData(dict=inputs['basis'])
        inputs['structure'] = self.ctx.structure_initial_primitive
        inputs['parameters'] = ParameterData(dict=inputs['parameters'])
        inputs['settings'] = ParameterData(dict=inputs['settings'])
        inputs['clean_workdir'] = Bool(False)
        inputs['max_iterations'] = Int(20)
        
        running = submit(SiestaBaseWorkChain, **inputs)
        self.report('launched SiestaBaseWorkChain<{}> in relaxation mode'.format(running.pid))
        
        return ToContext(workchain_relax=running)


    def run_bands(self):
        """
        Run the SiestaBaseWorkChain in scf+bands mode on the primitive cell of the relaxed input structure
        """
        self.report('Running bands calculation')

        try:
            structure = self.ctx.workchain_relax.out.output_structure
        except:
            self.abort_nowait('failed to get the output structure from the relaxation run')
            return
        
        self.ctx.structure_relaxed_primitive = structure


        inputs = dict(self.ctx.inputs)

        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(self.ctx.structure_relaxed_primitive)
        kpoints_mesh.set_kpoints_mesh_from_density(
            distance=self.ctx.protocol['kpoints_mesh_density'],
            offset=self.ctx.protocol['kpoints_mesh_offset'])

        bandskpoints = KpointsData()
        bandskpoints.set_cell(structure.cell, structure.pbc)
        bandskpoints.set_kpoints_path(kpoint_distance = 0.05)
        self.ctx.kpoints_path = bandskpoints

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        inputs['bandskpoints'] = self.ctx.kpoints_path           
        inputs['kpoints'] = kpoints_mesh
        inputs['structure'] = self.ctx.structure_relaxed_primitive
        inputs['parameters'] = ParameterData(dict=inputs['parameters'])
        inputs['basis'] = ParameterData(dict=inputs['basis'])
        inputs['settings'] = ParameterData(dict=inputs['settings'])
        
        running = submit(SiestaBaseWorkChain, **inputs)
        
        self.report('launching SiestaBaseWorkChain<{}> in scf+bands mode'.format(running.pid))
        
        return ToContext(workchain_bands=running)

    def run_results(self):
        """
        Attach the relevant output nodes from the band calculation to the workchain outputs
        for convenience
        """
        band_results = self.ctx.workchain_bands.out

        self.report('workchain succesfully completed'.format())
        self.out('scf_plus_band_parameters', band_results.output_parameters)
        self.out('bandstructure', band_results.bands_array)
        #self.out('remote_folder', calculation_band.out.remote_folder)
        #self.out('retrieved', calculation_band.out.retrieved)



