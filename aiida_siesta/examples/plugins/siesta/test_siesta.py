#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function

import sys
import os

from aiida.engine import submit
from aiida.orm import load_code
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.plugins import DataFactory

#  Siesta calculation on benzene molecule

PsfData = DataFactory('siesta.psf')
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
    "queue_name": "debug",
    "max_wallclock_seconds": 1700,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }
}

#TO DO:
# A Siesta executable compiled in serial mode might not work properly
# on a computer set up for MPI operation.
# This snippet can be used to check whether a code has been compiled
# with mpi support, and act accordingly
# For this to work, the user has to manually add the record in the
# database. In the verdi shell:
#
# code = load_node(code_PK)
# code.set_extra("mpi",True)
#code_mpi_enabled =  False
#try:
#    code_mpi_enabled =  code.get_extra("mpi")
#except AttributeError:
#    pass
#calc.set_withmpi(code_mpi_enabled)
#-----------------------------------------------------------------------

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
# Note the use of '.' in some entries. This will be fixed below.
# Note also that some entries have ':' as separator. This is not
# allowed in Siesta, and will be fixed by the plugin itself. The
# latter case is an unfortunate historical choice. It should not
# be used in modern scripts.
#
params_dict = {
    'xc-functional': 'LDA',
    'xc-authors': 'CA',
    'spin-polarized': True,
    'noncollinearspin': False,
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
pseudos_list = []
raw_pseudos = [("C.psf", ['C', 'Cred']), ("H.psf", 'H')]

for fname, kinds, in raw_pseudos:
    absname = os.path.realpath(
        os.path.join(os.path.dirname(__file__), "data/sample-psf-family",
                     fname))
    pseudo, created = PsfData.get_or_create(absname, use_first=True)
    if created:
        print("Created the pseudo for {}".format(kinds))
    else:
        print("Using the pseudo for {} from DB: {}".format(kinds, pseudo.pk))
    pseudos_list.append(pseudo)

#-----------------------------------------------------------------------

#
#--All the inputs of a Siesta calculations are listed in a dictionary--
#
inputs = {
    'structure': s,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'pseudos': {
        'C': pseudos_list[0],
        'Cred': pseudos_list[0],
        'H': pseudos_list[1],
    },
    'metadata': {
        'options': options,
        'label': "Benzene molecule",
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
