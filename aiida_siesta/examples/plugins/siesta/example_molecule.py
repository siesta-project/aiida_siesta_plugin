#!/usr/bin/env runaiida
import sys
import os.path as op

from aiida.engine import submit
from aiida.orm import load_code
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.plugins import DataFactory

#  Siesta calculation on benzene molecule

PsfData = DataFactory('pseudo.psf')
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
    codename = 'Siesta4.0.1@kelvin'

#
#------------------Code and computer options ---------------------------
#
code = load_code(codename)

options = {
#    "queue_name": "debug",
    "max_wallclock_seconds": 1700,
#    'withmpi': True,
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
alat = 15.  # angstrom
cell = [
    [
        alat,
        0.,
        0.,
    ],
    [
        0.,
        alat,
        0.,
    ],
    [
        0.,
        0.,
        alat,
    ],
]

# Note an atom tagged (for convenience) with a different label
s = StructureData(cell=cell)
s.append_atom(position=(0.000, 0.000, 0.468), symbols=['H'])
s.append_atom(position=(0.000, 0.000, 1.620), symbols=['C'])
s.append_atom(position=(0.000, -2.233, 1.754), symbols=['H'])
s.append_atom(position=(0.000, 2.233, 1.754), symbols=['H'])
s.append_atom(position=(0.000, -1.225, 2.327), symbols='C', name="Cred")
s.append_atom(position=(0.000, 1.225, 2.327), symbols=['C'])
s.append_atom(position=(0.000, -1.225, 3.737), symbols=['C'])
s.append_atom(position=(0.000, 1.225, 3.737), symbols=['C'])
s.append_atom(position=(0.000, -2.233, 4.311), symbols=['H'])
s.append_atom(position=(0.000, 2.233, 4.311), symbols=['H'])
s.append_atom(position=(0.000, 0.000, 4.442), symbols=['C'])
s.append_atom(position=(0.000, 0.000, 5.604), symbols=['H'])

elements = list(s.get_symbols_set())
#-----------------------------------------------------------------------

#
# ----------------------Parameters -------------------------------------
#
# Note the use of '.' is not allowed.
# Note also that ':' as separator is not allowed in Siesta.
# '-' is the suggested choice.
params_dict = {
    'xc-functional': 'LDA',
    'xc-authors': 'CA',
    'mesh-cutoff': '200.000 Ry',
    'max-scfiterations': 1000,
    'dm-numberpulay': 5,
    'dm-mixingweight': 0.050,
    'dm-tolerance': 1.e-4,
    'dm-mixscf1': True,
    'solution-method': 'diagon',
    'electronic-temperature': '100.000 K',
    'writeforces': True,
}

parameters = Dict(dict=params_dict)
#------------------------------------------------------------------------

#
# ---------------------Basis Set Info -----------------------------------
# The basis dictionary follows the 'parameters' convention
#
basis_dict = {
'pao-basistype':
'split',
'pao-splitnorm':
 0.150,
'pao-energyshift':
'0.020 Ry',
'%block pao-basis-sizes':
    """
C    SZP
Cred SZ
H    SZP
%endblock pao-basis-sizes""",
}

basis = Dict(dict=basis_dict)
#------------------------------------------------------------------------

#--------------------- Pseudopotentials ---------------------------------
#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
#
pseudos_dict = {}
raw_pseudos = [("C.psf", ['C', 'Cred']), ("H.psf", ['H'])]
for fname, kinds in raw_pseudos:
    absname = op.realpath(op.join(op.dirname(__file__), "../../fixtures/sample_psf", fname))
    pseudo = PsfData.get_or_create(absname)
    if not pseudo.is_stored:
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
    'pseudos' : pseudos_dict,
    'metadata': {
        'options': options,
        'label': "Benzene molecule",
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
