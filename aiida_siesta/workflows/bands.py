# -*- coding: utf-8 -*-
from __future__ import absolute_import
import six

from aiida.orm import Code
from aiida.orm import (Int, Str, Bool, Dict, StructureData, KpointsData,
                       BandsData)
from aiida.engine.launch import submit, run
from aiida.engine import WorkChain, ToContext, workfunction
from aiida.common.links import LinkType

from aiida_siesta.data.psf import get_pseudos_from_structure
##from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida_siesta.workflows.base import SiestaBaseWorkChain


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
            cls.run_bands,  # We can run this directly, a combined scf+bands
            cls.run_results,
        )

        # These will be inherited from the Base workchain outputs
        spec.output('bands', valid_type=BandsData)
        spec.output('output_structure',
                    valid_type=StructureData,
                    required=False)

        # This could be 'output_parameters', also inherited from the Base workchain, representing
        # the summary of energies, etc. coming out of the last SiestaCalculation run.
        # but its relevance is limited for a "band-structure" workflow. Left here for now for debugging.
        spec.output('output_parameters', valid_type=Dict)

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
            'code':
            self.inputs.code,
            'parameters': {},
            'settings': {},
            'options':
            Dict(dict={
                'resources': {
                    'num_machines': 1
                },
                'max_wallclock_seconds': 1800,
            }),
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
        inputs['clean_workdir'] = Bool(False)
        inputs['max_iterations'] = Int(20)

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report(
            'launched SiestaBaseWorkChain<{}> in relaxation mode'.format(
                running.pk))

        return ToContext(workchain_relax=running)

    def run_bands(self):
        """
        Run the SiestaBaseWorkChain in scf+bands mode on the primitive cell of the relaxed input structure
        """

        try:
            structure = self.ctx.workchain_relax.outputs.output_structure
        except:
            return self.exit_codes.ERROR_RELAXED_STRUCTURE_NOT_AVAILABLE

        # Do we need further refinement by Seekpath on this=? (eventually)?
        self.ctx.structure_relaxed_primitive = structure

        inputs = dict(self.ctx.inputs)

        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(
            self.ctx.structure_relaxed_primitive)
        kpoints_mesh.set_kpoints_mesh_from_density(
            distance=self.ctx.protocol['kpoints_mesh_density'],
            offset=self.ctx.protocol['kpoints_mesh_offset'])

        # For the band-structure kath, it is advised to use the
        # 'seekpath' method, but we try the 'legacy' for now.  In some
        # cases we might not want seekpath to change our structure.
        # Further support for this in the input to the workflow might
        # be needed.  (NOTE: If we ever optimize this workflow to
        # re-use the DM or H resulting from the execution of the Base
        # workflow, a change in structure would cause errors.)

        from aiida.tools import get_explicit_kpoints_path

        legacy_kpath_parameters = Dict(
            dict={
                'kpoint_distance':
                0.05  # In units of b1, b2, b3 (Around 20 points per side...)
            })
        seekpath_kpath_parameters = Dict(dict={
            'reference_distance': 0.02,
            'symprec': 0.0001
        })
        kpath_parameters = legacy_kpath_parameters

        result = get_explicit_kpoints_path(
            self.ctx.structure_relaxed_primitive,
            method='legacy',
            **kpath_parameters.get_dict())
        bandskpoints = result['explicit_kpoints']
        # The 'legacy' method presumably does not change the structure
        ## structure = result['primitive_structure']

        self.ctx.kpoints_path = bandskpoints

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        inputs['bandskpoints'] = self.ctx.kpoints_path
        inputs['kpoints'] = kpoints_mesh
        inputs['structure'] = self.ctx.structure_relaxed_primitive
        inputs['parameters'] = Dict(dict=inputs['parameters'])
        inputs['basis'] = Dict(dict=inputs['basis'])
        inputs['settings'] = Dict(dict=inputs['settings'])

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report(
            'launched SiestaBaseWorkChain<{}> in scf+bands mode'.format(
                running.pk))

        return ToContext(workchain_base_bands=running)

    def run_results(self):
        """
        Attach the relevant output nodes from the band calculation to the workchain outputs
        for convenience
        """

        from aiida.engine import ExitCode

        base_bands = self.ctx.workchain_base_bands

        #self.out('scf_plus_bands_summary', base_bands_results.output_parameters)
        #self.out('bands', base_bands_results.bands)
        for name, port in six.iteritems(self.spec().outputs):

            try:
                node = base_bands.get_outgoing(
                    link_label_filter=name).one().node
            except ValueError:
                if port.required:
                    self.report(
                        "the process spec specifies the output '{}' as required but was not an output of {}<{}>"
                        .format(name, base_bands, base_bands.pk))
            else:
                self.out(name, node)
                #self.report("attaching the node {}<{}> as '{}'"
                #            .format(node.__class__.__name__, node.pk, name))

        self.report('Bands workchain succesfully completed')
        return ExitCode(0)
