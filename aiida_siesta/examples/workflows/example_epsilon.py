#!/usr/bin/env runaiida

#Not required by AiiDA
import os.path as op
import sys

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import load_code
from aiida.orm import (Dict, StructureData, KpointsData)
from aiida_pseudo.data.pseudo.psf import PsfData
from aiida_siesta.workflows.epsilon import EpsilonWorkChain

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'

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
structure = StructureData(cell=cell)
structure.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
                      symbols=['Si'])
structure.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
                      symbols=['Si'])

#The parameters
parameters = {
        'xc-functional': 'LDA',
        'xc-authors': 'CA',
        'max-scfiterations': 50,
        'dm-numberpulay': 4,
        'dm-mixingweight': 0.3,
        'dm-tolerance': 1.e-3,
        'Solution-method': 'diagon',
        'electronic-temperature': '25 meV',
        'write-forces': True,
    }

relaxation = {
    'md-steps': 10
    }

#
# Use this for relaxation
#
parameters.update(relaxation)
parameters = Dict(dict=parameters)


#The basis set
basis = Dict(dict={
'pao-energy-shift': '300 meV',
'%block pao-basis-sizes': """
Si DZP
%endblock pao-basis-sizes""",
    })

# Optical calculation
optical = Dict(dict={
'optical-broaden': '0.5 eV',
'optical-energy-maximum': '20.0 eV',
'%block optical-mesh': """
 20 20 20
%endblock optical-mesh""",
'%block optical-vector': """
 1.0 0.0 0.0
%endblock optical-vector""",
    })

#The kpoints
kpoints = KpointsData()
kpoints.set_kpoints_mesh([4, 4, 4])

#The pseudopotentials
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
options = {
    "max_wallclock_seconds": 360,
    'withmpi': True,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 2,
    }
}

inputs = {
    'structure': structure,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'kpoints': kpoints,
    'pseudos': pseudos_dict,
    'options': Dict(dict=options),
    'optical': optical,
}

process = submit(EpsilonWorkChain, **inputs)
print("Submitted Epsilon workchain; ID={}".format(process.pk))
print(
    "For information about this workchain type: verdi process show {}".format(
        process.pk))
print("For a list of running processes type: verdi process list")
