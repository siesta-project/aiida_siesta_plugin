from __future__ import absolute_import
from aiida.engine import calcfunction, workfunction
from aiida.orm import Bool

@calcfunction
def get_protocol(struct, protocol, options, relax=Bool(False), autobands=Bool(False)):
    from aiida.orm import (KpointsData, Dict, StructureData)
    from aiida.tools import get_explicit_kpoints_path
    from aiida_siesta.data.psml import  get_pseudos_from_structure

    if protocol.value == 'standard':
        # Autobands
        structure=StructureData()
        if autobands.value:
            seekpath_param = dict={
                'reference_distance': 0.02,
                'symprec': 0.0001
            }
            result = get_explicit_kpoints_path(struct, **seekpath_param)
            structure = result['primitive_structure']
            bandskpoints = KpointsData()
            bandskpoints = result['explicit_kpoints']
            inputs={'bandskpoints' : bandskpoints}
        else:
            inputs={}
            structure = struct.clone()

        # Kpoints
        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(structure)
        kpoints_mesh.set_kpoints_mesh_from_density(distance=0.2,offset=[0., 0., 0.])

        #Parameters
        #First fixed
        parameters= {
            'max-scfiterations': 50,
            'dm-numberpulay': 4,
            'dm-mixingweight': 0.3,
            'dm-tolerance': 1.e-4,
            'md-max-force-tol': "0.02 eV/Ang",
            'solution-method': 'diagon',
            'electronic-temperature': '25 meV',
            'write-forces': True,
            'max-walltime': options["max_wallclock_seconds"]
        }
        #relaxations options
        if relax.value:
            parameters["md-type-of-run"]="cg"
            parameters["md-num-cg-steps"]=100
        #choice of meshcutoff including atomic_heuristics
        meshcutoff = 0.0
        min_meshcutoff = 100 # In Rydberg (!)
        atomic_heuristics = dict={
           'H': {'cutoff': 100 },
           'Si': {'cutoff': 100 }
                }
        for kind in structure.get_kind_names():
            try:
                cutoff = atom_heuristics[kind]['cutoff']
                meshcutoff = max(meshcutoff, cutoff)
            except:
                pass  # No problem. No heuristics, no info
        meshcutoff = max(min_meshcutoff,meshcutoff)  
        parameters["meshcutoff"]=meshcutoff

        #Basis
        basis={
           'pao-energy-shift': '100 meV',
           'pao-basis-size': 'DZP'
         }

        inputs.update({
                'structure' : structure,
                'basis' : Dict(dict=basis),
                'parameters' : Dict(dict=parameters),
                'kpoints' : kpoints_mesh,
                })

        return inputs

    else:
        return {}


@workfunction
def get_pseudo_p(structure, protocol):
    from aiida_siesta.data.psml import  get_pseudos_from_structure

    pseudos_dict = get_pseudos_from_structure(structure, 'nc-fr-04_pbe_stringent_psml')
    
    return pseudos_dict

