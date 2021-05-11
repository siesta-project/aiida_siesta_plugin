#!/usr/bin/env runaiida

#Not required by AiiDA
import os.path as op
import sys

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import load_code, load_node
from aiida.orm import (List, Dict, StructureData, KpointsData, Int, Float)
from aiida_siesta.data.psf import PsfData
from aiida_siesta.workflows.simplex_basis import SimplexBasisOptimization
from aiida_siesta.workflows.basis_optimization import BasisOptimizationWorkChain

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'

#The code
code = load_code(codename)

structure = load_node(62685)

#The parameters
parameters = Dict(
    dict={
        'meshcutoff': '500 Ry',
        'xc-functional': 'GGA',
        'xc-authors': 'PBE',
        'max-scfiterations': 4000,
        'scf-mixerhistory': 5,
        'scf-mixerweight': 0.1,
        'scf-dm-tolerance': 0.0001,
        'Solution-method': 'diagon',
        'electronic-temperature': '25 meV',
        'write-forces': True,
    })


#The kpoints
kpoints = KpointsData()
kpoints.set_kpoints_mesh([30, 30, 8])

#Resources
options = Dict(
    dict={
        "max_wallclock_seconds": 36000,
        "resources": {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        }
    })

#The submission
inputs = {
    'siesta_base': {
        'structure': structure,
        'parameters': parameters,
        'code': code,
        'kpoints': kpoints,
        'pseudos': {"C":load_node(23147)},
        'options': options
        },
    #'simplex': {
    #    },
    }

process = submit(BasisOptimizationWorkChain, **inputs)
#process = submit(TwoStepsBasisOpt, **inputs)
print(f"Submitted workchain; ID={process.pk}")
print(f"For information about this workchain type: verdi process show {process.pk}")
print("For a list of running processes type: verdi process list")
