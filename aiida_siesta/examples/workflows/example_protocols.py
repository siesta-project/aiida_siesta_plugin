"""
File showcasing the submission of a SiestaBaseWorkChain using the protocol system.
"""
import sys
from aiida_siesta.workflows.utils.base_inp_gen import BaseWorkchainInputsGenerator
from aiida.engine import submit
from aiida.orm import Dict, StructureData

try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'

#Three are the mandatory inputs: structure, a dictionary with the comp resources
#(calc_engines) and the protocol name

alat = 5.430  # angstrom
cell = [[0.5 * alat, 0.5 * alat, 0.,],
        [0., 0.5 * alat, 0.5 * alat,],
        [0.5 * alat, 0., 0.5 * alat,],
        ]
structure = StructureData(cell=cell)
structure.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
                      symbols=['Si'])
structure.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
                      symbols=['Si'])

calc_engines = {
     'siesta': {
         'code': codename, 
         'options': {
             'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1}, 
             "max_wallclock_seconds": 360, #'queue_name': 'DevQ', 'withmpi': True, 'account': "tcphy113c"
         }}}

protocol="standard_delta"

#Optionally, the bands calculation can be requested adding in inputs the `path_generator` parameter.
#Two options availables, the "legacy" option (no change of structure but wrong for some crystalin
#systems) and "seekpath"
path_generator = "seekpath"

#Optionally, the relaxation can be requested adding in inputs the `relaxation_type` parameter.
#At the moment the options "atoms_only", "variable_cell", "constant_volume" are available
relaxation_type = "atoms_only"


#We have now two options, either call the input generator:
inp_gen = BaseWorkchainInputsGenerator()
builder = inp_gen.get_builder(structure=structure, calc_engines=calc_engines, protocol=protocol)

#or access directly the method directly from the SiestaBaseWorkChain
builder = SiestaBaseWorkChain.get_filled_builder(structure=structure, calc_engines=calc_engines, protocol=protocol)

#The only advantage to use the input generator is that it has also some
#methods to guide the construction, like:
#inp_gen.get_relaxation_types()
#inp_gen.get_protocols_names()
#inp_gen.how_to_pass_computation_resources()

# As we get the builder suggested inputs, before submission any user has complete
# freedom to change something before submission.
# If no change is performed, just submitting the builder should still work and produce sensible results.
new_params = builder.parameters.get_dict()
new_params['max_scf_iterations'] = 200
builder.parameters = Dict(dict=new_params)

# Here we just submit the builder
process = submit(builder)
print("Submitted workchain; ID={}".format(process.pk))
print("For information about this workchain type: verdi process show {}".format(process.pk))
print("For a list of running processes type: verdi process list")

