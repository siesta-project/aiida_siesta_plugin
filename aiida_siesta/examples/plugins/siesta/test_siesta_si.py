#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import os.path as op
import sys

from aiida.engine import submit
from aiida.orm import load_code
from aiida.orm import (Dict, StructureData, KpointsData)
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.plugins import DataFactory

PsfData = DataFactory('siesta.psf')

try:
    dontsend = sys.argv[1]
    if dontsend == "--dont-send":
        submit_test = True
    elif dontsend == "--send":
        submit_test = False
    else:
        raise IndexError
except IndexError:
    print(("The first parameter can only be either "
           "--send or --dont-send"),
          file=sys.stderr)
    sys.exit(1)

try:
    codename = sys.argv[2]
except IndexError:
    codename = 'Siesta4.0.1@kelvin'

# Si diamond structure
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

# Si
# This was originally given in the "ScaledCartesian" format
#
structure = StructureData(cell=cell)
structure.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
                      symbols=['Si'])
structure.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
                      symbols=['Si'])

code = load_code(codename)

parameters = Dict(
    dict={
        'xc-functional': 'LDA',
        'xc-authors': 'CA',
        'max-scfiterations': 50,
        'dm-numberpulay': 4,
        'dm-mixingweight': 0.3,
        'dm-tolerance': 1.e-3,
        'Solution-method': 'diagon',
        'electronic-temperature': '25 meV',
        'md-typeofrun': 'cg',
        'md-numcgsteps': 3,
        'md-maxcgdispl': '0.1 Ang',
        'md-maxforcetol': '0.04 eV/Ang',
        'write-forces': True,
        # 'xml-write': True
    })

basis = Dict(
    dict={
        'pao-energy-shift': '300 meV',
        '%block pao-basis-sizes': """
Si DZP
%endblock pao-basis-sizes""",
    })

kpoints = KpointsData()
kpoints.set_kpoints_mesh([4, 4, 4])

pseudos_list = []
raw_pseudos = [("Si.psf", 'Si')]
for fname, kind in raw_pseudos:
    absname = op.realpath(
        op.join(op.dirname(__file__), "data/sample-psf-family", fname))
    pseudo, created = PsfData.get_or_create(absname, use_first=True)
    if created:
        print("\nCreated the pseudo for {}".format(kind))
    else:
        print("\nUsing the pseudo for {} from DB: {}".format(kind, pseudo.pk))
    pseudos_list.append(pseudo)

options = {
    "max_wallclock_seconds": 360,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }
}

inputs = {
    'structure': structure,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'kpoints': kpoints,
    'pseudos': {
        'Si': pseudos_list[0],
    },
    'metadata': {
        "label": "TestOnSiliconBulk",
        'options': options,
    }
}

if submit_test:
    inputs["metadata"]["dry_run"] = True
    inputs["metadata"]["store_provenance"] = False
    process = submit(SiestaCalculation, **inputs)
    #    subfolder, script_filename = calc.submit_test()
    print("Submited test for calculation (uuid='{}')".format(process.uuid))
    print("Check the folder submit_test for the result of the test")
# I could't find a way to access the actual folder (subfolder of submit_test)
# from the calculation node. So I can't print the exact location

else:
    process = submit(SiestaCalculation, **inputs)
    print("Submitted calculation; ID={}".format(process.pk))
    print("For information about this calculation type: verdi process show {}".
          format(process.pk))
    print("For a list of running processes type: verdi process list")
