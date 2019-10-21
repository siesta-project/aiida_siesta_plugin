#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import os.path as op
import sys
from aiida.tools import get_explicit_kpoints_path

#In this example we will calculate the band structure of Ge with SOC.
#Thanks to SeeK-path we can automatically generate the
#high symmetry points path where to calculate the bands.
#Alternatively, a manual list of k-points can be set.

################################################################

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

##-------------------Structure-----------------------------------
##Manually set the structure, all the quantities must be in Ang.
##Then, we pass through SeeK-path, to get the standardized cell,
##necessary for the automatic choice of the bands path.

alat = 5.65  # angstrom
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

s = StructureData(cell=cell)
s.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
              symbols=['Ge'])
s.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
              symbols=['Ge'])

elements = list(s.get_symbols_set())

seekpath_parameters = Dict(dict={
    'reference_distance': 0.02,
    'symprec': 0.0001
})
result = get_explicit_kpoints_path(s, **seekpath_parameters.get_dict())
structure = result['primitive_structure']

code = load_code(codename)

parameters = Dict(
    dict={
        'xc-functional': 'GGA',
        'xc-authors': 'PBE',
        'max-scfiterations': 50,
        'dm-numberpulay': 4,
        'Spin': "collinear",
        'dm-mixingweight': 0.3,
        'dm-tolerance': 1.e-3,
        'electronic-temperature': '25 meV'
    })

basis = Dict(
    dict={
        'pao-energy-shift': '300 meV',
        '%block pao-basis-sizes': """
Ge DZP
%endblock pao-basis-sizes""",
    })

kpoints = KpointsData()
kpoints.set_kpoints_mesh([4, 4, 4])

##-------------------K-points for bands --------------------
bandskpoints = KpointsData()

##.....Making use of SeeK-path for the automatic path......
##The choice of the distance between kpoints is in the call seekpath_parameters
##All high symmetry points included, labels already included
##This calls BandLine in siesta
bandskpoints = result['explicit_kpoints']

pseudos_list = []
raw_pseudos = [("Ge.psf", 'Ge')]
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
    "max_wallclock_seconds": 600,
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
    'bandskpoints': bandskpoints,
    'pseudos': {
        'Ge': pseudos_list[0],
    },
    'metadata': {
        "label": "Ge band-structure with collinear spin",
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

else:
    process = submit(SiestaCalculation, **inputs)
    print("Submitted calculation; ID={}".format(process.pk))
    print("For information about this calculation type: verdi process show {}".
          format(process.pk))
    print("For a list of running processes type: verdi process list")
