#!/usr/bin/env runaiida

#Not required by AiiDA
import os.path as op
import sys

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import load_code
from aiida.orm import (List, Dict, StructureData, KpointsData)
from aiida_siesta.data.psf import PsfData
from aiida_siesta.workflows.simplex_basis import SimplexBasisOptimization

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'

#The code
code = load_code(codename)

#The structure. Si diamond structure
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
        'xc-functional': 'LDA',
        'xc-authors': 'CA',
        'max-scfiterations': 4000,
        'dm-numberpulay': 4,
        'dm-mixingweight': 0.3,
        'dm-tolerance': 1.e-5,
        'Solution-method': 'diagon',
        'electronic-temperature': '25 meV',
        'write-forces': True,
    })

#The basis set
basis = Dict(
    dict={
        'pao-split-tail-norm':"T",
        '%block pao-basis': "\nSi   2\n n=3   0   2\n $0      0.0\n n=3   1   2\n $1      0.0\n%endblock pao-basis"
    })

#The kpoints
kpoints = KpointsData()
kpoints.set_kpoints_mesh([4, 4, 4])

#The pseudopotentials
pseudos_dict = {}
raw_pseudos = [("Si.psf", ['Si'])]
for fname, kinds in raw_pseudos:
    absname = op.realpath(
        op.join(op.dirname(__file__),
                "../plugins/siesta/data/sample-psf-family", fname))
    pseudo, created = PsfData.get_or_create(absname, use_first=True)
    if created:
        print("\nCreated the pseudo for {}".format(kinds))
    else:
        print("\nUsing the pseudo for {} from DB: {}".format(kinds, pseudo.pk))
    for j in kinds:
        pseudos_dict[j] = pseudo

#Resources
options = Dict(
    dict={
        "max_wallclock_seconds": 3600,
        "resources": {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        }
    })

#The submission
#All the inputs of a Siesta calculations are listed in a dictionary
#Note the different use of options compared to ../plugins/siesta/example_first.py
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
        'initial_variables': List(list=[4.0,6.0]),
        'upper_boundaries': List(list=[7.0,9.0]),
        'lower_boundaries': List(list=[2.0,3.0]),
       }
}

process = submit(SimplexBasisOptimization, **inputs)
print("Submitted workchain; ID={}".format(process.pk))
print(
    "For information about this workchain type: verdi process show {}".format(
        process.pk))
print("For a list of running processes type: verdi process list")
