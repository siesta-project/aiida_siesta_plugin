#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import sys
import os

__copyright__ = u"Copyright (c), This file is part of the AiiDA platform. For further information please visit http://www.aiida.net/. All rights reserved"
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file."
__version__ = "0.7.0"
__authors__ = "The AiiDA team."

#
# Prints benchmark info for a calculation
#

# Used to test the parent calculation
SiestaCalc = CalculationFactory('siesta.siesta') 

try:
    calc_id = sys.argv[1]
except IndexError:
    print(("Must provide as parameter the calc ID"), file=sys.stderr)
    sys.exit(1)


try:
    int(calc_id)
except ValueError:
    raise ValueError('Calc_id not an integer: {}'.format(calc_id))

calc = load_node(int(calc_id))

if isinstance(calc,SiestaCalc):

    print("Calculation status: '{}'".format(calc.get_state()), file=sys.stderr)
    print("Desc: {}".format(calc.description), file=sys.stderr)
    d=calc.out.output_parameters.get_dict()
    
    sys.stdout.write("{} ".format(int(calc_id)))
    try:
        no_u = d['no_u']
        nnz = d['nnz']
        print("No of orbitals: {} nnz: {}".format(no_u,nnz), file=sys.stderr)
        sys.stdout.write("{} {}".format(no_u,nnz))
    except:
        pass
    try:
        mesh = d['mesh']
        np = mesh[0]*mesh[1]*mesh[2]
        print("Mesh points: {} (total {})".format(mesh,np), file=sys.stderr)
        sys.stdout.write(" {} {}".format(mesh,np))
    except:
        pass
    try:
        nodes= d['siesta:Nodes']
        print("Number of nodes {}".format(nodes), file=sys.stderr)
        sys.stdout.write(" {} ".format(int(nodes)))
    except:
        pass
    try:
        code= str(calc.get_code()).split()[2]
        print("Code: {}".format(code), file=sys.stderr)
        sys.stdout.write(" {} ".format(code))
    except:
        pass
    try:
        computer= str(calc.get_computer()).split()[0]
        print("Computer: {}".format(computer), file=sys.stderr)
        sys.stdout.write(" {} ".format(computer))
    except:
        pass
    try:
        t_global = d['global_time']
        print("Total time: {}".format(t_global), file=sys.stderr)
    except:
        pass

    try:
        t_decomp = d['timing_decomposition']
        t_solver = t_decomp['compute_DM']
        t_H = t_decomp['setup_H']
        print("Total time, solver, setup_H: {} {} {}".format(t_global,t_solver,t_H), file=sys.stderr)
        sys.stdout.write(" {} {} {}\n".format(t_global,t_solver,t_H))
    except:
        pass
    
        sys.stdout.flush()

else:
    print(("Calculation should be a Siesta calculation."), file=sys.stderr)
    sys.exit(1)



