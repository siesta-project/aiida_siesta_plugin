#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-



import os.path as op
import sys
from aiida.tools import get_explicit_kpoints_path

#In this example we will calculate the band structure of Si.
#Thanks to SeeK-path we can automatically generate the
#high symmetry points path where to calculate the bands.
#Alternatively, a manual list of k-points can be set.

################################################################
import aiida
aiida.load_profile()

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

s = StructureData(cell=cell)
s.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
              symbols=['Si'])
s.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
              symbols=['Si'])

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

##-------------------K-points for bands --------------------
bandskpoints = KpointsData()
##Uncomment your favourite, two options:

##1)
##.....Making use of SeeK-path for the automatic path......
##The choice of the distance between kpoints is in the call seekpath_parameters
##All high symmetry points included, labels already included
##This calls BandLine in siesta
bandskpoints = result['explicit_kpoints']

##2)
##.....Only points, no labels.......
##Mandatory to set cell and pbc
##This calls BandsPoint
#kpp = [(0.500,  0.250, 0.750), (0.500,  0.500, 0.500), (0., 0., 0.)]
#bandskpoints.set_cell(structure.cell, structure.pbc)
#bandskpoints.set_kpoints(kpp)

#Note: The option to define a path touching specific kpoints
#for instance:
#kpp = [('W',  (0.500,  0.250, 0.750), 'L', (0.500,  0.500, 0.500), 40),
#        ('L', (0.500,  0.500, 0.500), 'G', (0., 0., 0.), 40)]
#Now is not easy to set. I'll study more on that

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
        'Si': pseudos_list[0],
    },
    'metadata': {
        "label": "TestOnSiliconBandsLines",
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
