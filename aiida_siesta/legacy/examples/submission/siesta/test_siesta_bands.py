#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import sys
import os

__copyright__ = u"Copyright (c), This file is part of the AiiDA platform. For further information please visit http://www.aiida.net/. All rights reserved"
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file."
__version__ = "0.7.0"
__authors__ = "The AiiDA team."

################################################################
UpfData = DataFactory('upf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')
RemoteData = DataFactory('remote')

# Used to test the parent calculation
SiestaCalc = CalculationFactory('siesta') 

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
    parent_id = sys.argv[2]
except IndexError:
    print >> sys.stderr, ("Must provide as second parameter the parent ID")
    sys.exit(1)


#####
# test parent

try:
    int(parent_id)
except ValueError:
    raise ValueError('Parent_id not an integer: {}'.format(parent_id))

parentcalc = Calculation.get_subclass_from_pk(parent_id)

#####

if isinstance(parentcalc,SiestaCalc):

    calc = parentcalc.create_restart()
    calc.label = "Test Siesta band structure"
    calc.description = "Test restart calculation with the Siesta code to get bands"

else:
    print >> sys.stderr, ("Parent calculation should be a Siesta "
                          "calculation.")
    sys.exit(1)

######

bandskpoints = KpointsData()
# Use the (possibly new) structure after a relaxation
bandskpoints.set_cell(calc.inp.structure.cell, calc.inp.structure.pbc)
bandskpoints.set_kpoints_path(kpoint_distance = 0.05)

calc.use_bandskpoints(bandskpoints)

#
# Make sure that we do a single-point calculation
# This is fragile, as it depends on these keys being
# in the same form as in the original calculation.

new_input_dict = calc.inp.parameters.get_dict()
new_input_dict['md-typeofrun'] = 'cg'
new_input_dict['md-numcgsteps'] = 0
new_input_dict['max-scfiterations'] = 0

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
    print "created calculation; calc=Calculation(uuid='{}') # ID: {}".format(
        calc.uuid,calc.dbnode.pk)
    calc.submit()
    print "submitted calculation; calc=Calculation(uuid='{}') # ID: {}".format(
        calc.uuid,calc.dbnode.pk)

