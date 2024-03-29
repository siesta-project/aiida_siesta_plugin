#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

#Not required by AiiDA
import os.path as op
import sys

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import Dict, KpointsData, StructureData, load_code
from aiida_pseudo.data.pseudo.psf import PsfData

from aiida_siesta.calculations.siesta import SiestaCalculation

##In alternative, Data and Calculation factories could be loaded.
##They containing all the data and calculation plugins:
#from aiida.plugins import DataFactory
#from aiida.plugins import CalculationFactory
#
#SiestaCalculation = CalculationFactory('siesta.siesta')
#PsfData = DataFactory('siesta.psf')
#StructureData = DataFactory('structure')
#...

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
    {
        'xc-functional': 'LDA',
        'xc-authors': 'CA',
        'max-scfiterations': 50,
        'dm-numberpulay': 4,
        'dm-mixingweight': 0.3,
        'dm-tolerance': 1.e-3,
        'Solution-method': 'diagon',
        'electronic-temperature': '25 meV',
        'write-forces': True,
    })

#The basis set
basis = Dict({
'pao-energy-shift': '300 meV',
'%block pao-basis-sizes': """
Si DZP
%endblock pao-basis-sizes""",
    })

#The kpoints
kpoints = KpointsData()
kpoints.set_kpoints_mesh([4, 4, 4])

#The pseudopotentials
pseudos_dict = {}
raw_pseudos = [("Si.psf", ['Si'])]
for fname, kinds in raw_pseudos:
    absname = op.realpath(op.join(op.dirname(__file__), "../../fixtures/sample_psf", fname))
    pseudo = PsfData.get_or_create(absname)
    if not pseudo.is_stored:
        print(f"\nCreated the pseudo for {kinds}")
    else:
        print(f"\nUsing the pseudo for {kinds} from DB: {pseudo.pk}")
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


#The submission
#All the inputs of a Siesta calculations are listed in a dictionary
inputs = {
    'structure': structure,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'kpoints': kpoints,
    'pseudos': pseudos_dict,
    'metadata': {
        "label": "TestOnSiliconBulk",
        'options': options,
    }
}

if submit_test:
    inputs["metadata"]["dry_run"] = True
    inputs["metadata"]["store_provenance"] = False
    process = submit(SiestaCalculation, **inputs)
    print(f"Submited test for calculation (uuid='{process.uuid}')")
    print("Check the folder submit_test for the result of the test")

else:
    process = submit(SiestaCalculation, **inputs)
    print(f"Submitted calculation; ID={process.pk}")
    print(f"For information about this calculation type: verdi process show {process.pk}")
    print("For a list of running processes type: verdi process list")

##An alternative is be to use the builder
#build=SiestaCalculation.get_builder()
#build.code=code
#build.structure=structure
#build.pseudos=pseudos_dict
#...
#build.metadata.options.resources = {'num_machines': 1 "num_mpiprocs_per_machine": 1}
#process = submit(builder)
