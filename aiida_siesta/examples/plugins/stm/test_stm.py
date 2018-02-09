#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
#
# Test for STM plugin
#
# Usage:
#
#    ./test_stm.py  {--send, --dont-send} codename remotedata_pk [height]
#
#  The codename argument is mandatory (the code is typically 'plstm', conforming to
#  the 'siesta.stm' plugin spec)
# 
#  The remotedata_pk argument is mandatory, and holds the PK of the object
#  representing the remote folder in which a LDOS file has been generated
#  (typically by a "relax+ldos" run of the SiestaBaseWorkchain workflow)
#
#  The height argument is optional. It is the height (in the z coordinate of
#  the LDOS box) at which the "image" is going to be computed.
#
#
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

#–-------------
try:
    codename = sys.argv[2]
except IndexError:
    print >> sys.stderr, ("Need a second parameter for the code name")
    sys.exit(1)

code = test_and_get_code(codename, expected_code_type='siesta.stm')
#–-------------
try:
    remotedata_pk = int(sys.argv[3])
except IndexError:
    print >> sys.stderr, ("Need a third parameter for the remotedata pk (LDOS)")
    sys.exit(1)
#–-------------
try:
    height = float(sys.argv[4])
except IndexError:
    height = 7.5
    
from aiida.orm.data.remote import RemoteData
remotedata = load_node(remotedata_pk)
#
#  Set up calculation object first
#
calc = code.new_calc()
calc.label = "Test STM"
calc.description = "STM calculation test"
#
#---- Attach remote data folder containing the LDOS
#
calc.use_parent_folder(remotedata)
#
#----Settings  -----------------------------
#
settings_dict={'additional_retrieve_list': []}
settings = ParameterData(dict=settings_dict)
calc.use_settings(settings)
#---------------------------------------------------

# Parameters ---------------------------------------------------
params_dict= {
    'z': height     # In Angstrom
}
parameters = ParameterData(dict=params_dict)
calc.use_parameters(parameters)
#
calc.set_max_wallclock_seconds(30*60) # 30 min
calc.set_resources({"num_machines": 1, "num_mpiprocs_per_machine": 1})


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

