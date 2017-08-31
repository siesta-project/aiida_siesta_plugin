#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia"

import sys
import os

from aiida.common.example_helpers import test_and_get_code
from aiida.common.exceptions import NotExistent

################################################################

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
    codename = sys.argv[2]
except IndexError:
    codename = 'plstm-4.0@rinaldo'

code = test_and_get_code(codename, expected_code_type='siesta.stm')
#
#  Set up calculation object first
#
calc = code.new_calc()
calc.label = "Test STM plugin. Benzene molecule"
calc.description = "Test calculation with the plstm code. Benzene molecule"

#
#----Settings first  -----------------------------
#
#settings_dict={'additional_retrieve_list': ['aiida.BONDS', 'aiida.EIG']}
#settings = ParameterData(dict=settings_dict)
#calc.use_settings(settings)
#---------------------------------------------------

#
# Parameters ---------------------------------------------------
#
params_dict= {
    'z': 3.0/0.529
}
parameters = ParameterData(dict=params_dict)
calc.use_parameters(parameters)
#
calc.set_max_wallclock_seconds(30*60) # 30 min

calc.set_resources({"num_machines": 1, "num_mpiprocs_per_machine": 1})
code_mpi_enabled =  False
try:
    code_mpi_enabled =  code.get_extra("mpi")
except AttributeError:
    pass
calc.set_withmpi(code_mpi_enabled)

from aiida.orm.data.remote import RemoteData
remotedata = load_node(2838)
calc.use_parent_folder(remotedata)

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

