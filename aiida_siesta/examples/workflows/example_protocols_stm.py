#!/usr/bin/env runaiida
"""
File showcasing the submission of a SiestaBaseWorkChain using the protocol system.
"""
import sys
from aiida_siesta.workflows.stm import SiestaSTMWorkChain
from aiida.engine import submit
from aiida.orm import Dict, StructureData, Float

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'
try:
    stmcodename = sys.argv[1]
except IndexError:
    stmcodename = 'STMhere@localhost'



# Structure -----------------------------------------
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

# Benzene molecule
# Note an atom tagged (for convenience) with a different label
s = StructureData(cell=cell)
s.append_atom(position=(0.468, 0.000, 0.000), symbols=['H'])
s.append_atom(position=(1.620, 0.000, 0.000), symbols=['C'])
s.append_atom(position=(1.754, -2.233, 0.000), symbols=['H'])
s.append_atom(position=(1.754, 2.233, 0.000), symbols=['H'])
s.append_atom(position=(2.327, -1.225, 0.000), symbols='C', name="Cred")
s.append_atom(position=(2.327, 1.225, 0.000), symbols=['C'])
s.append_atom(position=(3.737, -1.225, 0.000), symbols=['C'])
s.append_atom(position=(3.737, 1.225, 0.000), symbols=['C'])
s.append_atom(position=(4.311, -2.233, 0.000), symbols=['H'])
s.append_atom(position=(4.311, 2.233, 0.000), symbols=['H'])
s.append_atom(position=(4.442, 0.000, 0.000), symbols=['C'])
s.append_atom(position=(5.604, 0.000, 0.000), symbols=['H'])


calc_engines = {
     'siesta': {
         'code': codename, 
         'options': {'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1}, "max_wallclock_seconds": 3600 }
         },
     'stm': {
         'code': stmcodename,
         'options': {'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1}, "max_wallclock_seconds": 1360 }
         }
     }

protocol="standard_psml"

stm_mode = "constant-height"

stm_value = 1.6

inp_gen = SiestaSTMWorkChain.inputs_generator()
builder = inp_gen.get_filled_builder(s, calc_engines, protocol, stm_mode, stm_value)
#builder = inp_gen.get_builder(s, calc_engines, protocol, path_generator, relaxation_type, spin)

#The inputs generator (inp_gen) has also some
#methods to guide the construction, like:
print(inp_gen.get_relaxation_types())
print(inp_gen.get_protocol_names())
print(inp_gen.how_to_pass_computation_options())

# As we get the builder and not stored nodes, before submission any user has complete
# freedom to change something before submission.
# For instance the energy range for the construction of the LDOS.
builder.emin = Float(-6.4) #eV respect to Fermi energy

# Here we just submit the builder
process = submit(builder)
print("Submitted workchain; ID={}".format(process.pk))
print("For information about this workchain type: verdi process show {}".format(process.pk))
print("For a list of running processes type: verdi process list")

