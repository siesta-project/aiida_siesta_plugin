from __future__ import absolute_import
from aiida.engine import calcfunction, workfunction
from aiida.orm import Bool

#I would make it general, taking as input the 
#task! SiestaInputsGenerator("bands"), SiestaInputsGenerator("relax")
###
class SiestaRelaxationInputsGenerator:
    #Maybe we could create a data type, like for 
    #the pseudos and upload in the database.
    #This would be a SuperDataClass made of other data
    #classes (basis param pseudo)
    protocol_registy = dict={"standard" : "blabla",
                            "stringent" : "blublu"
                            }
    calc_types = dict={
            "relaxation": { 'code_plugin': 'siesta.siesta',
                            'description': 'These are calculations used for' 
                            'the main run of the code, computing the relaxation'
                          }
            }
    relax_types = dict={
            "atoms_only" : 
                "the latice shape and volume is fixed, only the athomic positions are relaxed",
             "variable_cell" : 
                "the lattice is relaxed together with the atomic coordinates. It allows"
                "to target hydro-static pressures or arbitrary stress tensors.",
            "constant_volume" : 
                "the cell volume is kept constant in a variable-cell relaxation: only" 
                "the cell shape andthe atomic coordinates are allowed to change.  Note that"
                "it does not make much sense tospecify a target stress or pressure in this" 
                "case, except for anisotropic (traceless) stresses"
              }


    @classmethod
    def get_calc_types(cls):
        return cls.calc_types.keys()

    @classmethod
    def get_calc_type_schema(cls, key):
        if key in cls.calc_types:
            return cls.calc_types[key]
        else:
            raise ValueError("wrong value of  calc_types")

    @classmethod
    def get_protocol_names(cls):
        return cls.protocol_registy.keys()

    @classmethod
    def get_default_protocol_name(cls):
        return "standard"

    @classmethod
    def get_protocol_info(cls, key):
        if key in cls.protocol_registy:
            return cls.protocol_registy[key]
        else:
            raise ValueError("wrong value of protocol")
    
    @classmethod
    def get_relaxation_types(cls):
        return cls.relax_types.keys()

    def get_builder(cls, structure, calc_engines, protocol,
            relaxation_type, threshold_forces=None, threshold_stress=None):

        from aiida_siesta.workflows.base import SiestaBaseWorkChain
        from aiida.orm import (Str, KpointsData, Dict, StructureData)
        from aiida.orm import load_code

        if protocol == 'stringent':

            kpoints_mesh = KpointsData()
            kpoints_mesh.set_cell_from_structure(structure)
            kpoints_mesh.set_kpoints_mesh_from_density(distance=0.2,offset=[0., 0., 0.])

            if not threshold_forces:
                threshold_forces=0.04

            if not threshold_stress:
                threshold_stress=0.006241509125883258 #1 GPa
            
            parameters= {
                'max-scfiterations': 50,
                'dm-numberpulay': 4,
                'dm-mixingweight': 0.3,
                'dm-tolerance': 1.e-4,
                'solution-method': 'diagon',
                'electronic-temperature': '25 meV',
                'write-forces': True,
                'max-walltime': calc_engines["relaxation"]["options"]["max_wallclock_seconds"],
                "md-type-of-run": "cg",
                "md-num-cg-steps": 100,
                "md-max-force-tol" : str(threshold_forces)+" eV/Ang",
                "md-max-stress-tol" : str(threshold_stress)+" eV/Ang**3"
                }
            #choice of meshcutoff including atomic_heuristics
            meshcutoff = 0.0
            min_meshcutoff = 100 # In Rydberg (!)
            atomic_heuristics = dict={
                'H': {'cutoff': 100 },
                'Si': {'cutoff': 101 }
                }
            for kind in structure.get_kind_names():
                try:
                    cutoff = atomic_heuristics[kind]['cutoff']
                    meshcutoff = max(meshcutoff, cutoff)
                except:
                    pass  # No problem. No heuristics, no info
            meshcutoff = max(min_meshcutoff,meshcutoff)  
            parameters["meshcutoff"]=str(meshcutoff)+" Ry"
            #relaxation type
            if  relaxation_type == "variable_cell":
                parameters["md-variable-cell"]= True
            if  relaxation_type == "constant_volume":
                parameters["md-constant-volume"]= True

            #Basis
            basis={
                'pao-energy-shift': '100 meV',
                'pao-basis-size': 'DZP'
            }

            #Pseudo fam
            pseudo_fam = "nc-fr-04_pbe_stringent_psml"

        else:
            if protocol != 'standard':
                import warnings
                warnings.warn("no protocol implemented with name {}, using default standard".format(protocol))
            
            kpoints_mesh = KpointsData()
            kpoints_mesh.set_cell_from_structure(structure)
            kpoints_mesh.set_kpoints_mesh_from_density(distance=0.3,offset=[0., 0., 0.])

            if not threshold_forces:
                threshold_forces=0.06

            if not threshold_stress:
                threshold_stress=0.01 #1.6 GPa
            
            parameters= {
                'max-scfiterations': 50,
                'dm-numberpulay': 4,
                'dm-mixingweight': 0.3,
                'dm-tolerance': 1.e-3,
                'solution-method': 'diagon',
                'electronic-temperature': '25 meV',
                'write-forces': True,
                'max-walltime': calc_engines["relaxation"]["options"]["max_wallclock_seconds"],
                "md-type-of-run": "cg",
                "md-num-cg-steps": 100,
                "md-max-force-tol" : str(threshold_forces)+" eV/Ang",
                "md-max-stress-tol" : str(threshold_stress)+" eV/Ang**3"
                }
            #choice of meshcutoff including atomic_heuristics
            meshcutoff = 0.0
            min_meshcutoff = 80 # In Rydberg (!)
            atomic_heuristics = dict={
                'H': {'cutoff': 100 },
                'Si': {'cutoff': 101 }
                }
            for kind in structure.get_kind_names():
                try:
                    cutoff = atomic_heuristics[kind]['cutoff']
                    meshcutoff = max(meshcutoff, cutoff)
                except:
                    pass  # No problem. No heuristics, no info
            meshcutoff = max(min_meshcutoff,meshcutoff)  
            parameters["meshcutoff"]=str(meshcutoff)+" Ry"
            #relaxation type
            if  relaxation_type == "variable_cell":
                parameters["md-variable-cell"]= True
            if  relaxation_type == "constant_volume":
                parameters["md-constant-volume"]= True

            #Basis
            basis={
                'pao-energy-shift': '100 meV',
                'pao-basis-size': 'DZ'
            }
            
            #Pseudo fam
            pseudo_fam = "sample_psf_family"

        builder = SiestaBaseWorkChain.get_builder()
        builder.structure = structure
        builder.basis = Dict(dict=basis)
        builder.parameters = Dict(dict=parameters)
        builder.kpoints = kpoints_mesh
        builder.pseudo_family = Str(pseudo_fam)
        builder.options = Dict(dict=calc_engines["relaxation"]["options"])
        builder.code = code = load_code(calc_engines["relaxation"]["code"])

        return builder

