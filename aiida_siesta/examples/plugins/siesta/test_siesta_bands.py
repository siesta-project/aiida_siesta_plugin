#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import sys
import os

__copyright__ = u"Copyright (c), This file is part of the AiiDA platform. For further information please visit http://www.aiida.net/. All rights reserved"
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file."
__version__ = "0.7.0"
__authors__ = "The AiiDA team."


# This script will get as an input the PK of a previous calculation 
# and it will restart it to get the band structure.

# ISSUE: to generate standard k-point path for the band structure,
#        SeeK-path uses a standardized cell. In normal situations the 
#        correct approach to use this tool is the following:
#        1) You first find the standardized primitive cell with SeeK-path 
#           (returned in output) together with the k-point coordinates and 
#           suggested band path
#        2) You then run all your calculations using the SeeK-path
#           standardized primitive cell
#        For example see test_siesta_si_bands.py or test_siesta_cif_bands.py
#        Here the situation is more complicated as we start from a previous
#        calculation, so an already fixed structure. At the moment I don't
#        use SeeK-path, manual kpoints for band.

# Usage:
#         ./test_siesta_bands.py {--send, --dont-send} PK_of_original_calculation


################################################################
UpfData = DataFactory('upf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')
RemoteData = DataFactory('remote')

# Used to test the parent calculation
SiestaCalc = CalculationFactory('siesta.siesta') 

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
## Use the (possibly new) structure after a relaxation

bandskpoints = KpointsData()

kpp = [('W',  (0.500,  0.250, 0.750), 'L', (0.500,  0.500, 0.500), 40),
        ('L', (0.500,  0.500, 0.500), 'G', (0., 0., 0.), 40)]
from aiida.tools.data.array.kpoints import legacy
fs=legacy.get_explicit_kpoints_path(kpp)
bandskpoints.set_cell(calc.inp.structure.cell, calc.inp.structure.pbc)
bandskpoints.set_kpoints(fs[3])
bandskpoints.labels=fs[4]


## IN CASE WE WANT AUTOMATIC K POINT PATH.
## SEEKPATH USES A STANDARD STRUCTURE, CAN WE CHANGE THE STRUCTURE IN THE RESTART??????
## IF WE COULD DO THAT, THEN:
#from aiida.tools import get_explicit_kpoints_path
#seekpath_parameters = ParameterData(dict={'reference_distance': 0.02,'symprec': 0.0001})
#result=get_explicit_kpoints_path(calc.inp.structure, **seekpath_parameters.get_dict())
#newstructure = result['primitive_structure']
#set newstructure as structure of the calcultion!! How??
#bandskpoints=result['explicit_kpoints']

calc.use_bandskpoints(bandskpoints)

#
# Make sure that we do a single-point calculation
# This is fragile, as it depends on these keys being
# in the same form as in the original calculation.

new_input_dict = calc.inp.parameters.get_dict()
new_input_dict['md-typeofrun'] = 'cg'
new_input_dict['md-numcgsteps'] = 0
new_input_dict['max-scfiterations'] = 50

calc.use_parameters(ParameterData(dict=new_input_dict))

if submit_test:
    subfolder, script_filename = calc.submit_test()
    print "Test_submit for calculation (uuid='{}')".format(calc.uuid)
    print "Submit file in {}".format(os.path.join(os.path.relpath(subfolder.abspath),script_filename))
else:
    calc.store_all()
    print "created calculation; calc=Calculation(uuid='{}') # ID: {}".format(calc.uuid,calc.dbnode.pk)
    calc.submit()
    print "submitted calculation; calc=Calculation(uuid='{}') # ID: {}".format(calc.uuid,calc.dbnode.pk)
