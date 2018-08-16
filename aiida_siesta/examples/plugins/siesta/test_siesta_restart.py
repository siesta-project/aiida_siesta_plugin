#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

__copyright__ = u"Copyright (c), 2017, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"

# This script will restart a calculation that ended in a FAILED state due
# to lack of scf convergence in the allotted number of iterations.
# ISSUE: is it possible to increase the number of iteration from here???
# Usage:
#         ./test_siesta_restart.py {--send, --dont-send} PK_of_failed_calculation

import sys, os
from aiida.orm import load_node

ParameterData = DataFactory('parameter')

try:
    dontsend = sys.argv[1]
    if dontsend == "--dont-send":
        submit_test = True
    elif dontsend == "--send":
        submit_test = False
    else:
        raise IndexError
except IndexError:
    print >> sys.stderr, ("The first parameter can only be either "
                          "--send or --dont-send")
    sys.exit(1)

try:
    PK = int(sys.argv[2])
except IndexError:
    print >> sys.stderr, ("The second parameter must be the PK of a calculation")
    sys.exit(1)

    
c = load_node(PK)
print "Restarting calculation (uuid='{}')".format(c.uuid)
print "Calculation status: '{}'".format(c.get_state())
print " "
calc = c.create_restart(force_restart=True)

new_input_dict = c.inp.parameters.get_dict()
new_input_dict['max-scfiterations'] = 50
calc.use_parameters(ParameterData(dict=new_input_dict))

if submit_test:
    subfolder, script_filename = calc.submit_test()
    print "Test_submit for calculation (uuid='{}')".format(
        calc.uuid)
    print "Submit file in {}".format(os.path.join(
        os.path.relpath(subfolder.abspath),
        script_filename
        ))
else:
    calc.store_all()
    print "created calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.dbnode.pk)
    calc.submit()
    print "submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.dbnode.pk)

