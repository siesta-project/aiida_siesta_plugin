#!/usr/bin/env runaiida

#Not required by AiiDA
import os.path as op
import sys

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import load_code
from aiida.orm import (Dict, StructureData, KpointsData)
from aiida.orm import SinglefileData

from aiida_siesta.data.psf import PsfData
from aiida_siesta.workflows.exchange_barrier import ExchangeBarrierWorkChain

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'

code = load_code(codename)

from aiida.orm import StructureData
cell = [[ 9.9678210000 ,  0.0000000000 ,    0.0000000000,],
        [-4.9839090000 ,  8.6323860000 ,    0.0000000000,],
        [ 0.0000000000 ,  0.0000000000 ,  35.0000000000,],
]
s = StructureData(cell=cell)
s.append_atom(position=( 1.6613036650053632    , 0.959153999040846 ,17.5 ),symbols='Mg',name="Mg") #1
s.append_atom(position=( 3.333333331578814e-07 ,  1.918307998081692,17.5 ),symbols='O',name="O")   #2
s.append_atom(position=( 4.983911666666668     , 6.714078001918309 ,17.5 ),symbols='Mg',name="Mg") #3
s.append_atom(position=( 3.322608334994638     , 7.673232000959154 ,17.5 ),symbols='O',name="O")   #4
s.append_atom(position=( 1.6613046600214534    ,  6.714078001918309,17.5 ),symbols='Mg',name="Mg") #5
s.append_atom(position=( 1.3283494233462534e-06, 7.673232000959154 ,17.5 ),symbols='O',name="O")   #6
s.append_atom(position=(-1.6613023366559396    ,  6.714078001918309,17.5 ),symbols='Mg',name="Mg") #7
s.append_atom(position=(-3.32260566832797      , 7.673232000959154 ,17.5 ),symbols='O',name="O")   #8
s.append_atom(position=( 6.645214669989273     , 3.836615996163384 ,17.5 ),symbols='Mg',name="Mg") #9
s.append_atom(position=( 4.983911338317244     , 4.79576999520423  ,17.5 ),symbols='O',name="O")   #10
s.append_atom(position=( 3.3226076633440593    , 3.836615996163384 ,17.5 ),symbols='Mg',name="Mg") #11
s.append_atom(position=( 1.6613043316720297    , 4.79576999520423  ,17.5 ),symbols='O',name="O")   #12
s.append_atom(position=( 6.666666663157628e-07 , 3.836615996163384 ,17.5 ),symbols='Mg',name="Mg") #13
s.append_atom(position=(-1.6613026650053635    , 4.79576999520423  ,17.5 ),symbols='O',name="O")   #14
s.append_atom(position=( 8.30651766832797      , 0.959153999040846 ,17.5 ),symbols='Mg',name="Mg") #15
s.append_atom(position=( 6.6452143366559415    , 1.918307998081692 ,17.5 ),symbols='O',name="O_a")   #16
s.append_atom(position=( 4.983910661682756     , 0.959153999040846 ,17.5 ),symbols='Mg',name="Mg") #17
s.append_atom(position=( 3.3226073300107264    , 1.918307998081692 ,17.5 ),symbols='O',name="O_b")  #18

structure = s

# Exchange O_a and O_b above (note python base-0 convention for indexes)
i1 = 15
i2 = 17
migration_direction = [ 0.0, 0.0, 1.0 ]    # Z direction

# Lua script
absname = op.abspath(
        op.join(op.dirname(__file__), "../plugins/siesta/lua_scripts/neb.lua"))
n_images_in_script=5
lua_script = SinglefileData(absname)


# Parameters: very coarse for speed of test

parameters = dict={
   "mesh-cutoff": "50 Ry",
   "dm-tolerance": "0.001",
   "DM-NumberPulay ":  "3",
   "DM-History-Depth":  "0",
   "SCF-Mixer-weight":  "0.02",
   "SCF-Mix":   "density",
   "SCF-Mixer-kick":  "35",
   "MD-VariableCell":  "F",
   "MD-MaxCGDispl":  "0.3 Bohr",
   "MD-MaxForceTol":  " 0.04000 eV/Ang"
    }

# All other atoms are fixed (...)
constraints = dict={
    "%block Geometry-Constraints":
    """
    atom [ 1 -- 15 ]
    %endblock Geometry-Constraints"""
    }

#
# Use this for constraints
#
parameters.update(constraints)
#
neb_parameters = Dict(dict=parameters)

# Extra parameter for end-point relaxation
relaxation = dict={
    'md-steps': 10
    }

parameters.update(relaxation)
endpoint_parameters = Dict(dict=parameters)

    
# The basis set
basis = Dict(dict={
'pao-energy-shift': '300 meV',
'%block pao-basis-sizes': """
Mg SZ
O SZ
O_a SZ
O_b SZ
%endblock pao-basis-sizes""",
    })


# The kpoints
kpoints_endpoints = KpointsData()
kpoints_endpoints.set_kpoints_mesh([1,1,1])

kpoints_neb = KpointsData()
kpoints_neb.set_kpoints_mesh([1,1,1])

# The pseudopotentials
pseudos_dict = {}
raw_pseudos = [("Mg.psf", ['Mg']), ("O.psf", ['O','O_a','O_b'])]
for fname, kinds in raw_pseudos:
    absname = op.realpath(
        op.join(op.dirname(__file__), "../plugins/siesta/data/sample-psf-family", fname))
    pseudo, created = PsfData.get_or_create(absname, use_first=True)
    if created:
        print("\nCreated the pseudo for {}".format(kinds))
    else:
        print("\nUsing the pseudo for {} from DB: {}".format(kinds, pseudo.pk))
    for j in kinds:
        pseudos_dict[j]=pseudo

# Resources
options = {
    "max_wallclock_seconds": 3600,
    'withmpi': True,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 2,
    }
}

#
# For finer-grained compatibility with script
options_neb = {
    "max_wallclock_seconds": 7200,
    'withmpi': True,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 2,
    }
}

endpoint_inputs= {
    'parameters': endpoint_parameters,
    'code': code,
    'basis': basis,
    'kpoints': kpoints_endpoints,
    'pseudos': pseudos_dict,
    'options': Dict(dict=options)
}


inputs = {

    'initial_structure': structure,
    'first_index':      Int(i1),
    'second_index':     Int(i2),
    'migration_direction': List(list=migration_direction),
    'n_images': Int(n_images_in_script),
    
    'initial': endpoint_inputs,
    'final': endpoint_inputs,

    'neb': {
        'neb_script': lua_script,
        'parameters': neb_parameters,
        'code': code,
        'basis': basis,
        'kpoints': kpoints_neb,
        'pseudos': pseudos_dict,
        'options': Dict(dict=options_neb)
    },
        
}

process = submit(ExchangeBarrierWorkChain, **inputs)
print("Submitted ExchangeBarrier workchain; ID={}".format(process.pk))
print("For information about this workchain type: verdi process show {}".format(process.pk))
print("For a list of running processes type: verdi process list")
