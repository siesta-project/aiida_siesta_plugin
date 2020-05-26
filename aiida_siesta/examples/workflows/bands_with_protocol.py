"""
File showcasing the submission of a SiestaBaseWorkChain in order to calculate bands dispersions using
inputs generated through the method `get_builder` of the class `SiestaBandsInputsGenerator`.
"""
import sys
from aiida_siesta.workflows.utils.bandsinputs import SiestaBandsInputsGenerator
from aiida.engine import submit
from aiida.orm import Dict, StructureData

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'

calc_engines = {
     'bands': {
         'code': codename, 
         'options': {
             'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1}, 
             "max_wallclock_seconds": 360, #'queue_name': 'DevQ', 'withmpi': True, 'account': "tcphy113c"
         }}}

protocol="standard_delta"

#The authomatic generation of the high symmetry path (defining where to calculate the bands)
#is usually performed by seekpath. Seekpath requires to change the input structure to follow
#some convension. If the user doesn't want the cell to be touched, the "legacy" option
#is available.
path_generator = "seekpath"

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
structure = StructureData(cell=cell)
structure.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
                      symbols=['Si'])
structure.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
                      symbols=['Si'])

#We now create an instance of the class, but it is not
#neccesarry as get_builder is classmethod
bands_inp_gen = SiestaBandsInputsGenerator()

# This is the main call: we get a builder for the `SiestaBaseWorkChain`, pre-filled, 
# with unstored nodes (unless they are taken from the DB, e.g. pseudos)
builder = bands_inp_gen.get_builder(
    structure=structure, calc_engines=calc_engines, protocol=protocol, 
    path_generator=path_generator)

# The user now received a builder with suggested inputs, but before submission he/she has complete
# freedom to change any of them.
# NOTE: The changes are code-specific and optional.
# If no change is performed, just submitting the builder should still work and produce sensible results.
new_params = builder.parameters.get_dict()
new_params['max_scf_iterations'] = 200
builder.parameters = Dict(dict=new_params)

# Here we just submit the builder
process = submit(builder)
print("Submitted workchain; ID={}".format(process.pk))
print("For information about this workchain type: verdi process show {}".format(process.pk))
print("For a list of running processes type: verdi process list")

