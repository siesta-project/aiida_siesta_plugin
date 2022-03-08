#!/usr/bin/env runaiida

# LUA PATH MUST BE PASSED AS SECOND OPTION!!!!!!

# An example driver for the Siesta NEB Base workchain, using
# also the feature to include ghost sites (floating orbitals)
# in a transparent manner.
# The ghost information is specified in a special entry
# 'floating_sites' in the basis dictionary.
# The image data in the ../fixtures/neb-data concerns
# only the physical atoms, but the workchain takes care to add
# the floating sites when it generates the NEB folder-data for the
# Siesta calculation. (Note that this feature is not available when
# a NEB calculation is setup directly in a Siesta calculation).

# In output a node `neb_output_package` is produced. Use the
# function plot_neb in aiida_siesta.utils.neb to plot the results.

#Not required by AiiDA
import os.path as op
import sys

#AiiDA classes and functions
from aiida.engine import submit
from aiida.orm import load_code
from aiida.orm import (Dict, StructureData, KpointsData)
from aiida.orm import TrajectoryData, SinglefileData
from aiida_pseudo.data.pseudo.psf import PsfData
from aiida_siesta.utils.xyz_utils import get_structure_list_from_folder
from aiida_siesta.workflows.neb_base import SiestaBaseNEBWorkChain


try:
    codename = sys.argv[1]
    load_code(codename)
except (IndexError, NotExistent):
    print(("The first parameter must be the code to use. Hint: `verdi code list`."),file=sys.stderr)
    sys.exit(1)

try:
    lua_elements_path = sys.argv[2]
except IndexError:
    print(("The second parameter must be the path to the lua scripts in the flos library."),file=sys.stderr)
    print(("Look at the docs for more info. Library can be found at https://github.com/siesta-project/flos"),file=sys.stderr)
    sys.exit(1)


#The code
code = load_code(codename)

cell = [[15.0, 00.0 , 00.0,],
        [00.0, 15.0 , 00.0,],
        [00.0, 00.0 , 15.0,],
        ]
s = StructureData(cell=cell)
s.append_atom(position=( 0.000,  0.000,  0.000 ),symbols=['O']) #1
s.append_atom(position=( 0.757,  0.586,  0.000 ),symbols=['H']) #2
s.append_atom(position=(-0.757,  0.586,  0.000 ),symbols=['H']) #3 
s.append_atom(position=( 0.000,  3.500,  0.000 ),symbols=['O']) #4
s.append_atom(position=( 0.757,  2.914,  0.000 ),symbols=['H']) #5
s.append_atom(position=(-0.757,  2.914,  0.000 ),symbols=['H']) #6

image_structure_list = get_structure_list_from_folder("../fixtures/neb_data", s)
_kinds_raw = [ k.get_raw() for k in image_structure_list[0].kinds ]

path_object = TrajectoryData(image_structure_list)
path_object.set_attribute('kinds', _kinds_raw)

# Lua script
absname = op.abspath(op.join(op.dirname(__file__), "../fixtures/lua_scripts/neb.lua"))
lua_script = SinglefileData(absname)


#The parameters
#
# NOTE that we put "by hand" an extra constraint
# on the ghost atom. Without it, the NEB algorithm
# would likely not converge, as the magnitude of the forces on
# ghosts bear no relation to the rest...
#
parameters = Dict(dict={
   "mesh-cutoff": "50 Ry",
   "dm-tolerance": "0.0001",
   "DM-NumberPulay ":  "3",
   "DM-History-Depth":  "0",
   "SCF-Mixer-weight":  "0.02",
   "SCF-Mix":   "density",
   "SCF-Mixer-kick":  "35",
   "MD-VariableCell":  "F",
   "MD-MaxCGDispl":  "0.3 Bohr",
   "MD-MaxForceTol":  " 0.04000 eV/Ang",
    "%block Geometry-Constraints":
    """
    atom [1 -- 4]
    atom 7
    %endblock Geometry-Constraints"""
    })

# Basis set info, including floating sites. Note that
# O_top will get a default DZP basis set. If needed,
# it can be specified in the block.
basis = Dict(dict={
  'floating_sites': [ {"name":'O_top', "symbols":'O', "position":(-0.757,  0.586,  2.00 ) } ],
  '%block PAO-Basis':
    """
 O                     2                    # Species label, number of l-shells
 n=2   0   2                         # n, l, Nzeta
   3.305      2.510
   1.000      1.000
 n=2   1   2 P   1                   # n, l, Nzeta, Polarization, NzetaPol
   3.937      2.542
   1.000      1.000
H                     1                    # Species label, number of l-shells
 n=1   0   2 P   1                   # n, l, Nzeta, Polarization, NzetaPol
   4.828      3.855
   1.000      1.000

    %endblock PAO-Basis""",
})


# The pseudopotentials
# Lacking an automatic procedure to recognize the 'O' symbol in 'O_top'
# (considered unsafe), we need to add explicitly the 'O_top' kind name
# to the 'O.psf' association list.
# Still to check: use of pseudopotential families
pseudos_dict = {}
raw_pseudos = [ ("H.psf", ['H']),("O.psf", ['O', 'O_top'])]
for fname, kinds in raw_pseudos:
    absname = op.realpath(op.join(op.dirname(__file__), "../fixtures/sample_psf", fname))
    pseudo = PsfData.get_or_create(absname)
    if not pseudo.is_stored:
        print("\nCreated the pseudo for {}".format(kinds))
    else:
        print("\nUsing the pseudo for {} from DB: {}".format(kinds, pseudo.pk))
    for j in kinds:
        pseudos_dict[j]=pseudo

# Resources and other options
options = Dict(dict={
    "max_wallclock_seconds": 3600,
    'withmpi': True,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 2,
    },
    "environment_variables":{"LUA_PATH":lua_elements_path},
})


# The workchain submission

inputs = {
    'starting_path': path_object,
    'neb_script': lua_script,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'pseudos': pseudos_dict,
    'options': options
}

process = submit(SiestaBaseNEBWorkChain, **inputs)
print("Submitted Siesta NEB Base workchain; ID={}".format(process.pk))
print("For information type: verdi process show {}".format(process.pk))
print("For a list of running processes type: verdi process list")

