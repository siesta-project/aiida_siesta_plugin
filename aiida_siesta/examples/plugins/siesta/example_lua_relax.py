#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

#LUA PATH MUST BE PASSED AS THIRD OPTION!!!!!!

import os.path as op
import sys

from aiida.engine import submit
from aiida.orm import Group, SinglefileData, load_code
from aiida.plugins import DataFactory

from aiida_siesta.calculations.siesta import SiestaCalculation

Dict = DataFactory('core.dict')
KpointsData = DataFactory('core.array.kpoints')
StructureData = DataFactory('core.structure')

try:
    dontsend = sys.argv[1]
    if dontsend == "--dont-send":
        submit_test = True
    elif dontsend == "--send":
        submit_test = False
    else:
        raise IndexError
except IndexError:
    print(("The first parameter can only be either --send or --dont-send."),file=sys.stderr)
    sys.exit(1)

try:
    codename = sys.argv[2]
    load_code(codename)
except (IndexError, NotExistent):
    print(("The second parameter must be the code to use. Hint: `verdi code list`."),file=sys.stderr)
    sys.exit(1)

try:
    lua_elements_path = sys.argv[3]
except IndexError:
    print(("The third parameter must be the path to the lua scripts in the flos library."),file=sys.stderr)
    print(("Look at the docs for more info. Library can be found at https://github.com/siesta-project/flos"),file=sys.stderr)
    sys.exit(1)
    #lua_elements_path = "/home/ebosoni/flos/?.lua;/home/ebosoni/flos/?/init.lua;"


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
    },
    "environment_variables":{"LUA_PATH":lua_elements_path},
}
#
#-------------------------- Settings ---------------------------------
#
settings_dict = {'additional_retrieve_list': ['aiida.BONDS', 'aiida.EIG']}
settings = Dict(settings_dict)
#
# Structure -----------------------------------------
#
alat = 10.0  # angstrom
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

# Water molecule
# One of the H atoms is sligthy moved

s = StructureData(cell=cell)
s.append_atom(position=(0.000, 0.000, 0.00), symbols=['O'])
s.append_atom(position=(0.757, 0.586, 0.00), symbols=['H'])
s.append_atom(position=(-0.780, 0.600, 0.00), symbols=['H'])

# ----------------------Parameters -------------------------------------

params_dict = {
    'xc-functional': 'LDA',
    'xc-authors': 'CA',
    'mesh-cutoff': '100.000 Ry',
    'max-scfiterations': 30,
    'dm-numberpulay': 4,
    'dm-mixingweight': 0.1,
    'dm-tolerance': 1.e-4,
    'md-maxcgdispl': '0.200 bohr',
    'md-maxforcetol': '0.020 eV/Ang',
}

parameters = Dict(params_dict)
#------------------------------------------------------------------------
#
# No basis set spec in this calculation (default)
#
#--------------------- Pseudopotentials ---------------------------------
#
# FIXME: The family name is hardwired
#
family = Group.get(label='psf_family')
pseudos_dict = family.get_pseudos(structure=s)
#-----------------------------------------------------------------------
# Lua script for relaxation
#
absname = op.abspath(op.join(op.dirname(__file__), "../../fixtures/lua_scripts/relax_geometry_lbfgs.lua"))
lua_script = SinglefileData(absname)
#
#
inputs = {
    'lua': { 'script': lua_script},
    'structure': s,
    'parameters': parameters,
    'code': code,
    'pseudos': pseudos_dict,
    'metadata': {
        'options': options,
        'label': "Water molecule relaxation with LBFGS"
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
