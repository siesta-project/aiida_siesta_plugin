#!/usr/bin/env runaiida

#Not required by AiiDA
import os.path as op
import sys

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import load_code, load_node
from aiida.orm import (Str, List, Dict, StructureData, KpointsData, Int, Float)
from aiida_pseudo.data.pseudo.psf import PsfData
from aiida_siesta.workflows.simplex_basis import SimplexBasisOptimization
from aiida_siesta.workflows.two_steps_optimization import TwoStepsBasisOpt

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
        'meshcutoff': '100 Ry',
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

#The basis set 'pao-split-tail-norm':"T",
basis = Dict(
    dict={
        '%block pao-basis': "\nSi   2\n n=3   0   2\n 4.99376      $sz2 \n n=3   1   2 P 1\n 6.2538      $pz2 \n%endblock pao-basis"
    })

#The kpoints
kpoints = KpointsData()
kpoints.set_kpoints_mesh([8, 8, 8])

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
        'basis': basis,
        'kpoints': kpoints,
        'pseudos': pseudos_dict,
        'options': options
        },
    'simplex': {
     #   'max_iters': Int(4),
        'output_name': Str("basis_enthalpy"),
        'variables_dict': Dict(dict={
            "sz2":[2.0,4.8,3.0],
            "pz2":[2.0,6.0,3.0]
            }),
        },
    #'macrostep':{'lambda_scaling_factor': Float(0.2)},
    }

process = submit(SimplexBasisOptimization, **inputs)
#process = submit(TwoStepsBasisOpt, **inputs)
print(f"Submitted workchain; ID={process.pk}")
print(f"For information about this workchain type: verdi process show {process.pk}")
print("For a list of running processes type: verdi process list")
