#!/usr/bin/env runaiida

#Not required by AiiDA
import os.path as op
import sys
import io

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import load_code
from aiida.orm import (Int, Float, Dict, StructureData, KpointsData)
from aiida_pseudo.data.pseudo.psf import PsfData
from aiida_siesta.workflows.eos import EqOfStateFixedCellShape

# This example shows the use of the EqOfStateFixedCellShape
# Requires a working aiida profile and the set up of
# a code (it submits the WorkChain to the daemon).
# To run it: runaiida example_eos.py codename
# The inputs are the same of ../plugins/siesta/example_first.py
# with the addition of an optional "volume_per_atom",
# the starting volume around which the EoS is calculated

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'Siesta-4.0.2@kay'

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
        'max-scfiterations': 40,
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
        'pao-energy-shift': '300 meV',
        '%block pao-basis-sizes': """
Si DZP
%endblock pao-basis-sizes""",
    })

#The kpoints
kpoints = KpointsData()
kpoints.set_kpoints_mesh([14, 14, 14])

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
options = Dict(
    dict={
        "max_wallclock_seconds": 360,
        #'withmpi': True,
        #'account': "tcphy113c",
        #'queue_name': "DevQ",
        "resources": {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        }
    })

#The submission
#All the inputs of a Siesta calculations are listed in a dictionary
#Note the different use of options compared to ../plugins/siesta/example_first.py
inputs = {
    'structure': structure,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'kpoints': kpoints,
    'pseudos': pseudos_dict,
    'options': options,
    'volume_per_atom': Float(19),
    'batch_size': Int(7) #selects the number of volumes to run at the same time
}

process = submit(EqOfStateFixedCellShape, **inputs)
print("Submitted workchain; ID={}".format(process.pk))
print(
    "For information about this workchain type: verdi process show {}".format(
        process.pk))
print("For a list of running processes type: verdi process list")
