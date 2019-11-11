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

# Calculation on Iron, collinear spin polarization applied

################################################################

PsfData = DataFactory('siesta.psf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')

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
    codename = 'siesta4.0.1@parsons'

code = test_and_get_code(codename, expected_code_type='siesta.siesta')
#
#  Set up calculation object first
#
calc = code.new_calc()
calc.label = "Siesta-Fe-bulk-spin"
calc.description = "Test Siesta calculation. Fe bulk spin-polarized"

#
# Structure -----------------------------------------
#
# BCC
alat = 2.87 # angstrom
cell = [[alat/2, alat/2, alat/2,],
        [alat/2, -alat/2, alat/2,],
        [alat/2, alat/2, -alat/2,],
       ]

#
s = StructureData(cell=cell)
s.append_atom(position=(0.000,0.000,0.000),symbols=['Fe'])

elements = list(s.get_symbols_set())
calc.use_structure(s)
#-------------------------------------------------------------
#
# Parameters ---------------------------------------------------
#
# Note the use of '.' in some entries. This will be fixed below.
# Note also that some entries have ':' as separator. This is not
# allowed in Siesta, and will be fixed by the plugin itself. The
# latter case is an unfortunate historical choice. It should not
# be used in modern scripts.
#
params_dict= {
'xc.functional': 'GGA',
'xc.authors': 'PBE',
'spin-polarized': True,
'noncollinearspin': False,
'mesh-cutoff': '150.000 Ry',
'max-scfiterations': 40,
'dm-numberpulay': 4,
'dm-mixingweight': 0.1,
'dm-tolerance': 1.e-3,
'electronic-temperature': '300.000 K'
}
#
# Sanitize, as '.' is not kosher for the database handlers
#
params_dict = { k.replace('.','-') :v for k,v in params_dict.iteritems() }
#
parameters = ParameterData(dict=params_dict)
calc.use_parameters(parameters)
#
kpoints = KpointsData()
kpoints_mesh = 6
kpoints.set_kpoints_mesh([kpoints_mesh,kpoints_mesh,kpoints_mesh])
calc.use_kpoints(kpoints)
#
# Pseudopotentials ----------------------------------------------
#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
# Families support is not yet available for this.
#
raw_pseudos = [ ("Fe.psf", 'Fe')]

for fname, kinds, in raw_pseudos:
    absname = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                            "data",fname))
    pseudo, created = PsfData.get_or_create(absname,use_first=True)
    if created:
        print "Created the pseudo for {}".format(kinds)
    else:
        print "Using the pseudo for {} from DB: {}".format(kinds,pseudo.pk)
        
    # Attach pseudo node to the calculation
    calc.use_pseudo(pseudo,kind=kinds)
#-------------------------------------------------------------------

calc.set_max_wallclock_seconds(30*60) # 30 min

calc.set_resources({"num_machines": 1})
#------------------

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

