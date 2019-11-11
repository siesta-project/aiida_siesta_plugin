#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function

# See LICENSE and Contributors

# This script will restart a calculation that ended without
# scf convergence in the allotted number of iterations, or without geometry convergence.
#
# Usage:
#         ./test_siesta_restart.py {--send, --dont-send} PK_of_failed_calculation

import sys, os

from aiida.engine import submit
from aiida.orm import load_code
from aiida.common import NotExistent
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.plugins import DataFactory

Dict = DataFactory('dict')

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
                          "--send or --dont-send"), file=sys.stderr)
    sys.exit(1)

try:
    PK = int(sys.argv[2])
except IndexError:
    print(("The second parameter must be the PK of a calculation"), file=sys.stderr)
    sys.exit(1)

    
c = load_node(PK)

print("Restarting calculation (uuid='{}')".format(c.uuid))
print("Is excepted?: '{}'".format(c.is_excepted))
#
# Note that this use of 'failed' does not correspond to the process state.
#
print("Is failed?: '{}'".format(c.is_failed))
if c.is_failed:
   print("Exit code: '{}'".format(c.exit_status))
   print("'{}'".format(c.exit_message))
print(" ")

restart=c.get_builder_restart()

newpar=restart.parameters.clone()
newpar.attributes["max-scf-iterations"]= 50
restart.parameters=newpar

if c.outputs.output_parameters.attributes["variable_geometry"]:
   restart.structure=c.outputs.output_structure

# The most important line. The presence of
# parent_calc_folder triggers the 'restart' operations
# in the plugin, such as  the copy of the .DM and the
# addition of use-save-dm to the parameters

restart.parent_calc_folder=c.outputs.remote_folder

if submit_test:

    m=restart["metadata"]

    # An attempt to clone the 'metadata' section of the object raised this error:
    # Error: AttributeError: 'ProcessBuilderNamespace' object has no attribute 'clone'
    # It seems that one can simply set this attributes by hand:

    m.dry_run = True
    m.store_provenance = False
    print("Dry run without storing provenance. See the submit_test directory")

else:
    print("Calculation restarted")


submit(restart)

