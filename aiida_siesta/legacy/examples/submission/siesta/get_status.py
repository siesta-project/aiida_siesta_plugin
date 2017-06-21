#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import sys
import os

__copyright__ = u"Copyright (c), This file is part of the AiiDA platform. For further information please visit http://www.aiida.net/. All rights reserved"
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file."
__version__ = "0.7.0"
__authors__ = "The AiiDA team."

#
# Prints warnings for a calculation
#

# Used to test the parent calculation
SiestaCalc = CalculationFactory('siesta.siesta') 

try:
    calc_id = sys.argv[1]
except IndexError:
    print >> sys.stderr, ("Must provide as parameter the calc ID")
    sys.exit(1)


try:
    int(calc_id)
except ValueError:
    raise ValueError('Calc_id not an integer: {}'.format(calc_id))

#calc = Calculation.get_subclass_from_pk(calc_id)
calc = load_node(int(calc_id))
#####

if isinstance(calc,SiestaCalc):

    print "Calculation status: '{}'".format(calc.get_state())
##??     print "Desc: {}".format(calc.description)
    d=calc.out.output_parameters.get_dict()
    
    try:
        print "Warnings: {}".format(d['warnings'])
    except:
        pass

else:
    print >> sys.stderr, ("Calculation should be a Siesta calculation.")
    sys.exit(1)



