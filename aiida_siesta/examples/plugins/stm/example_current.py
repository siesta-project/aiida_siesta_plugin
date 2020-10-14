#!/usr/bin/env runaiida
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

import sys

from aiida.engine import submit
from aiida.orm import load_code
from aiida.orm import load_node
from aiida.orm import Dict, Str, Float
from aiida_siesta.calculations.stm import STMCalculation
from aiida.plugins import DataFactory

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
#
#------------------Code and computer options ---------------------------
#
try:
    codename = sys.argv[2]
except IndexError:
    codename = 'Siesta4.0.1@kelvin'

code = load_code(codename)
#
#–-------------
try:
    remotedata_pk = int(sys.argv[3])
except IndexError:
    print(("Need a third parameter for the remotedata pk (LDOS)"),
          file=sys.stderr)
    sys.exit(1)

remotedata = load_node(remotedata_pk)

#–-------------
try:
    cur = float(sys.argv[4])
except IndexError:
    cur = 0.3

#
options = {
    #"queue_name": "debug",
    "max_wallclock_seconds": 1700,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }
}
#
# Parameters ---------------------------------------------------
#
settings_dict = {}
settings = Dict(dict=settings_dict)
#
#--All the inputs of a Siesta calculations are listed in a dictionary--
#
inputs = {
    'settings': settings,
    'value': Float(cur),
    'mode': Str("constant-current"),
    'code': code,
    'ldos_folder': remotedata,
    'metadata': {
        'options': options,
        'label': "STM test",
    }
}

if submit_test:
    inputs["metadata"]["dry_run"] = True
    inputs["metadata"]["store_provenance"] = False
    process = submit(STMCalculation, **inputs)
    print("Submited test for calculation (uuid='{}')".format(process.uuid))
    print("Check the folder submit_test for the result of the test")

else:
    process = submit(STMCalculation, **inputs)
    print("Submitted calculation; ID={}".format(process.pk))
    print("For information about this calculation type: verdi process show {}".
          format(process.pk))
    print("For a list of running processes type: verdi process list")
