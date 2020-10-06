#!/usr/bin/env runaiida

import os.path as op
import sys
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.engine import submit

# Script to restart a calculation, for instance the one obtained with
# example_scf_fail.py or example_geom_fail.py
# Usage:
#         ./example_restart.py {--send, --dont-send} PK_of_calculation_to_restart


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

try:
    pke = sys.argv[2]
    pk = int(pke)
except:
    print(("The second parameter is the pk of the calculation "
           "you want to restart, integer, mandatory!"),
          file=sys.stderr)
    sys.exit(1)

#Loading calculation we want to restart (old_cal),
#must be second parameter passed to the script
try:
    g=load_node(pk)
    if (g.process_class == SiestaCalculation):
        print(("Restarting calculation {}").format(pk))
    else:
        raise
except:
    print(("PK you passed is not a valid PK. "
           "Allowed are CalcJobNode with .process_class == SiestaCalculation"))
    sys.exit(1)


#Set up the a new calculation with all
#the inputs of the old one (we use the builder)
restart=g.get_builder_restart()

#The inputs of old_calc attched to restart are
#already stored!!! You can not modify them straight away
#If you want to change something you make a clone
#and reassign it to the builder.
#Here we change max-scfiterations for example
newpar=restart.parameters.clone()
newpar.attributes["max-scfiterations"]=32
##In case of geometry fail, change th number of md steps:
#newpar.attributes["md-numcgsteps"]=20
restart.parameters=newpar


#We need to take care here of passing the
#output geometry of old_calc to the new calculation
if g.outputs.output_parameters.attributes["variable_geometry"]:
   restart.structure=g.outputs.output_structure

#The most important line. The presence of
#parent_calc_folder triggers the real restart
#meaning the copy of the .DM and the
#addition of dm-use-saved-dm to the parameters
restart.parent_calc_folder=g.outputs.remote_folder

if submit_test:
    restart.metadata.dry_run=True
    restart.metadata.store_provenance=False
    process = submit(restart)
    print("Submited test for calculation (uuid='{}')".format(process.uuid))
    print("Check the folder submit_test for the result of the test")

else:
    process = submit(restart)
    print("Submitted calculation; ID={}".format(process.pk))
    print("For information about this calculation type: verdi process show {}".
          format(process.pk))
    print("For a list of running processes type: verdi process list")

