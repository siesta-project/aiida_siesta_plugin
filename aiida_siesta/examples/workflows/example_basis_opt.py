#!/usr/bin/env runaiida

#Not required by AiiDA
import os.path as op
import sys

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import load_code, load_node
from aiida.orm import (List, Dict, Bool, StructureData, KpointsData, Int, Float)
from aiida_pseudo.data.pseudo.psf import PsfData
from aiida_siesta.workflows.simplex_basis import SimplexBasisOptimization
from aiida_siesta.workflows.basis_optimization import BasisOptimizationWorkChain

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'

#The code
code = load_code(codename)

#Structure
alat = 5.430  # angstrom
cell = [
    [
        0.5 * alat,
        0.5 * alat,
        0.,
    ],
    [
        0.,
        0.5 * alat,
        0.5 * alat,
    ],
    [
        0.5 * alat,
        0.,
        0.5 * alat,
    ],
]
#The atom positions were originally given in the "ScaledCartesian" format
#but standard for aiida structures is Cartesian in Angstrom
structure = StructureData(cell=cell)
structure.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
                      symbols=['Si'])
structure.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
                      symbols=['Si'])

#The parameters
parameters = Dict(
    dict={
        'xc-functional': "GGA",
        'xc-authors': "PBE",
        'max-scf-iterations': 200,
        'scf-mixer-history': 5,
        'scf-mixer-weight': 0.1,
        'scf-dm-tolerance': 1.e-5,
        'solution-method': 'diagon',
        'electronic-temperature': '25 meV',
        'write-forces': True,
        'mesh-cutoff': '100 Ry',
        #'scf-dm-tolerance': 0.0001,
    })


#The kpoints
kpoints = KpointsData()
kpoints.set_kpoints_mesh([8, 8, 8])

#Resources
options = Dict(
    dict={
        "max_wallclock_seconds": 36000,
        "resources": {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        }
    })

#Pseudos
pseudos_dict = {}
raw_pseudos = [("Si.psf", ['Si'])]
for fname, kinds in raw_pseudos:
    absname = op.realpath(op.join(op.dirname(__file__), "../fixtures/sample_psf", fname))
    pseudo = PsfData.get_or_create(absname)
    if not pseudo.is_stored:
        print("\nCreated the pseudo for {}".format(kinds))
    else:
        print("\nUsing the pseudo for {} from DB: {}".format(kinds, pseudo.pk))
    for j in kinds:
        pseudos_dict[j]=pseudo


#The submission
inputs = {
    'siesta_base': {
        'structure': structure,
        'parameters': parameters,
        'code': code,
        'kpoints': kpoints,
        'pseudos': pseudos_dict,
        'options': options
        },
    'simplex': {"max_iters" : Int(10), "tolerance_function": Float(0.1)},
    'optimization_schema': {'global_split_norm':Bool(True)}
    }

process = submit(BasisOptimizationWorkChain, **inputs)
#process = submit(TwoStepsBasisOpt, **inputs)
print(f"Submitted workchain; ID={process.pk}")
print(f"For information about this workchain type: verdi process show {process.pk}")
print("For a list of running processes type: verdi process list")
