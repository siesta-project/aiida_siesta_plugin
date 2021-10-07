#!/usr/bin/env runaiida
import os.path as op
import sys

#In this example we will calculate the band structure of Si.
#Thanks to SeeK-path we can automatically generate the
#high symmetry points path where to calculate the bands.
################################################################

from aiida.engine import submit
from aiida.orm import load_code
from aiida.orm import (Dict, StructureData, KpointsData)
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida.orm import Group
from aiida.tools import get_explicit_kpoints_path

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


#The code
codename = 'siesta-school--MaX-1.3.0-1@localhost'
code = load_code(codename)

#Structure
#We pass through SeeK-path, to get the standardized cell,
#necessary for the automatic choice of the bands path.
alat = 5.430  # angstrom
cell = [[0.5 * alat, 0.5 * alat, 0.],
        [0., 0.5 * alat, 0.5 * alat],
        [0.5 * alat,0.,0.5 * alat]]
s = StructureData(cell=cell)
s.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
              symbols=['Si'])
s.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
              symbols=['Si'])
seekpath_parameters = Dict(dict={'reference_distance': 0.02})
result = get_explicit_kpoints_path(s, **seekpath_parameters.get_dict())
structure = result['primitive_structure']

#The parameters
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
        'write-forces': True,
        'mesh-cutoff': "200 Ry"
    })

#The basis
basis = Dict(dict={
'pao-energy-shift': '100 meV',
'%block pao-basis-sizes': """
Si DZP
%endblock pao-basis-sizes""",
    })

#The kpoints mesh
kpoints = KpointsData()
kpoints.set_kpoints_mesh([11, 11, 11])

##-------------------K-points for bands --------------------
bandskpoints = KpointsData()
bandskpoints = result['explicit_kpoints']

#The pseudopotentials
family = Group.get(label='PseudoDojo/0.4/PBE/SR/standard/psml')
pseudos_dict = family.get_pseudos(structure=structure)

#Resources
options = {
    "max_wallclock_seconds": 600,
#    "withmpi" : True,
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
    'pseudos': pseudos_dict,
    'metadata': {
        "label": "TestOnSiliconBands",
    }
}

#SiestaCalulation
inputs['metadata']['options'] = options

#SiestaBaseWorkChain
#inputs['options'] = Dict(dict=options)

if submit_test:
    inputs["metadata"]["dry_run"] = True
    inputs["metadata"]["store_provenance"] = False
    process = submit(SiestaCalculation, **inputs)
    print("Submited test for calculation (uuid='{}')".format(process.uuid))
    print("Check the folder submit_test for the result of the test")

else:
    process = submit(SiestaCalculation, **inputs)
    #process = submit(SiestaBaseWorkChain, **inputs)
    print("Submitted calculation; ID={}".format(process.pk))
    print("For information about this calculation type: verdi process show {}".
          format(process.pk))
    print("For a list of running processes type: verdi process list")

