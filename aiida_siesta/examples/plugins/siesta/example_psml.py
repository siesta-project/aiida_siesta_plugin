#!/usr/bin/env runaiida

import sys
import os.path as op

from aiida.engine import submit
from aiida.orm import load_code
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.plugins import DataFactory

# Version of siesta supporting psml files is needed to run 
# this example
##########################################################

PsmlData = DataFactory('siesta.psml')
Dict = DataFactory('dict')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')

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
#
try:
    codename = sys.argv[2]
except IndexError:
    codename = 'SiestaMax@kelvin'

#
#------------------Code and computer options ---------------------------
#
code = load_code(codename)

options = {
#    "queue_name": "debug",
    "max_wallclock_seconds": 1700,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }
}

#
#-------------------------- Settings ---------------------------------
#
settings_dict = {'additional_retrieve_list': ['aiida.BONDS', 'aiida.EIG']}
settings = Dict(dict=settings_dict)
#---------------------------------------------------------------------

#
#-------------------------- Structure --------------------------------
#
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
s.append_atom(
    position=(0.000 * alat, 0.000 * alat, 0.000 * alat), symbols=['Si'])
s.append_atom(
    position=(0.250 * alat, 0.250 * alat, 0.250 * alat), symbols=['Si'])

#-----------------------------------------------------------------------

#
# ----------------------Parameters -------------------------------------
#
# Note the use of '.' in some entries. This will be fixed below.
# Note also that some entries have ':' as separator. This is not
# allowed in Siesta, and will be fixed by the plugin itself. The
# latter case is an unfortunate historical choice. It should not
# be used in modern scripts.
#
params_dict = {
    'xc-functional': 'LDA',
    'xc-authors': 'CA',
    'spin': 'SO',
    'mesh-cutoff': '200.000 Ry',
    'max-scfiterations': 1000,
    'dm-numberpulay': 5,
    'dm-mixingweight': 0.050,
    'dm-tolerance': 1.e-4,
    'dm-mixscf1': True,
    'negl-nonoverlap-int': False,
    'solution-method': 'diagon',
    'electronic-temperature': '100.000 K',
    'md-typeofrun': 'cg',
    'md-numcgsteps': 2,
    'md-maxcgdispl': '0.200 bohr',
    'md-maxforcetol': '0.050 eV/Ang',
    'writeforces': True,
    'writecoorstep': True,
    'write-mulliken-pop': 1,
}

parameters = Dict(dict=params_dict)
#------------------------------------------------------------------------

#
# ---------------------Basis Set Info -----------------------------------
# The basis dictionary follows the 'parameters' convention
#
basis_dict = {
    'pao-basistype': 'split',
    'pao-splitnorm': 0.150,
    'pao-energyshift': '0.020 Ry'
}

basis = Dict(dict=basis_dict)
#------------------------------------------------------------------------

#--------------------- Pseudopotentials ---------------------------------
#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
#
pseudos_dict = {}
raw_pseudos = [("Si.psml", ['Si'])]
for fname, kinds in raw_pseudos:
    absname = op.realpath(
        op.join(op.dirname(__file__), "data/sample-psml-family", fname))
    pseudo, created = PsmlData.get_or_create(absname, use_first=True)
    if created:
        print("\nCreated the pseudo for {}".format(kinds))
    else:
        print("\nUsing the pseudo for {} from DB: {}".format(kinds, pseudo.pk))
    for j in kinds:
        pseudos_dict[j]=pseudo


#-----------------------------------------------------------------------

#
#--All the inputs of a Siesta calculations are listed in a dictionary--
#
inputs = {
    'structure': s,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'pseudos': pseudos_dict,
    'metadata': {
        'options': options,
        'label': "Si crystal with PSML pseudos",
    }
}

if submit_test:
    inputs["metadata"]["dry_run"] = True
    inputs["metadata"]["store_provenance"] = False
    process = submit(SiestaCalculation, **inputs)
    print("Submited test for calculation (uuid='{}')".format(process.uuid))
    print("Check the folder submit_test for the result of the test")

else:
    process = submit(SiestaCalculation, **inputs)
    print("Submitted calculation; ID={}".format(process.pk))
    print("For information about this calculation type: verdi process show {}".
          format(process.pk))
    print("For a list of running processes type: verdi process list")
