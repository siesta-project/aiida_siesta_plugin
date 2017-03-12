#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import sys
import os
import time
import subprocess

from aiida.common.example_helpers import test_and_get_code
from aiida.common.exceptions import NotExistent

expected_free_energy = -1015.377304
timeout_secs = 5*60.
queue = None

################################################################

PsfData = DataFactory('psf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')

codename = 'siesta@torquessh'

code = test_and_get_code(codename, expected_code_type='siesta')
#
#  Set up calculation object first
#
calc = code.new_calc()
calc.label = "Test Siesta. Benzene molecule"
calc.description = "Test calculation with the Siesta code. Benzene molecule"

#
#----Settings first  -----------------------------
#
settings_dict={'additional_retrieve_list': ['aiida.BONDS', 'aiida.EIG']}
settings = ParameterData(dict=settings_dict)
calc.use_settings(settings)
#---------------------------------------------------

#
# Structure -----------------------------------------
#
alat = 15. # angstrom
cell = [[alat, 0., 0.,],
        [0., alat, 0.,],
        [0., 0., alat,],
       ]

# Benzene molecule
# Note an atom tagged (for convenience) with a different label
#
s = StructureData(cell=cell)
s.append_atom(position=(0.000,0.000,0.468),symbols=['H'])
s.append_atom(position=(0.000,0.000,1.620),symbols=['C'])
s.append_atom(position=(0.000,-2.233,1.754),symbols=['H'])
s.append_atom(position=(0.000,2.233,1.754),symbols=['H'])
s.append_atom(position=(0.000,-1.225,2.327),symbols='C',name="Cred")
s.append_atom(position=(0.000,1.225,2.327),symbols=['C'])
s.append_atom(position=(0.000,-1.225,3.737),symbols=['C'])
s.append_atom(position=(0.000,1.225,3.737),symbols=['C'])
s.append_atom(position=(0.000,-2.233,4.311),symbols=['H'])
s.append_atom(position=(0.000,2.233,4.311),symbols=['H'])
s.append_atom(position=(0.000,0.000,4.442),symbols=['C'])
s.append_atom(position=(0.000,0.000,5.604),symbols=['H'])

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
'xc.functional': 'LDA',
'xc.authors': 'CA',
'spin:polarized': True,
'noncollinearspin': False,
'mesh-cutoff': '200.000 Ry',
'max-scfiterations': 1000,
'dm-numberpulay': 5,
'dm-mixingweight': 0.050,
'dm-tolerance': 1.e-4,
'dm-mixscf1': True,
'negl-nonoverlap-int': False,
'solution-method': 'diagon',
'electronic-temperature': '100.000 K',
'md-typeofrun': 'cg',
'md-numcgsteps': 2,
'md-maxcgdispl': '0.200 bohr',
'md-maxforcetol': '0.050 eV/Ang',
'writeforces': True,
'writecoorstep': True,
'xml-write': True,
'write-mulliken-pop': 1,
'%block example-block': """
first line
second line    """,
}
#
# Sanitize, as '.' is not kosher for the database handlers
#
params_dict = { k.replace('.','-') :v for k,v in params_dict.iteritems() }
#
parameters = ParameterData(dict=params_dict)
calc.use_parameters(parameters)
#
#----------------------------------------------------------
#
# Basis Set Info ------------------------------------------
# The basis dictionary follows the 'parameters' convention
#
basis_dict = {
'pao-basistype': 'split',
'pao-splitnorm': 0.150,
'pao-energyshift': '0.020 Ry',
'%block pao-basis-sizes' :"""
C    SZP
Cred SZ
H    SZP  """,
}
#
basis_dict = { k.replace('.','-') :v for k,v in  basis_dict.iteritems() }
#
basis = ParameterData(dict=basis_dict)
calc.use_basis(basis)
#--------------------------------------------------------------


# Pseudopotentials ----------------------------------------------
#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
# Families support is not yet available for this.
#
raw_pseudos = [ ("C.psf", ['C', 'Cred']),
                ("H.psf", 'H')]

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

## For remote codes, it is not necessary to manually set the computer,
## since it is set automatically by new_calc
#computer = code.get_remote_computer()
#calc = code.new_calc(computer=computer)

calc.set_max_wallclock_seconds(30*60) # 30 min

calc.set_resources({"num_machines": 1})
calc.set_withmpi(False)
#------------------
queue = None
# calc.set_queue_name(queue)
#------------------

#from aiida.orm.data.remote import RemoteData
#calc.set_outdir(remotedata)

calc.store_all()
print "created calculation; calc=Calculation(uuid='{}') # ID={}".format(
    calc.uuid,calc.dbnode.pk)
calc.submit()
print "submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
    calc.uuid,calc.dbnode.pk)

#####################################

print "Wating for end of execution..."
start_time = time.time()
exited_with_timeout = True
while time.time() - start_time < timeout_secs:
    time.sleep(15) # Wait a few seconds
    

    # print some debug info, both for debugging reasons and to avoid
    # that the test machine is shut down because there is no output

    print "#"*78
    print "####### TIME ELAPSED: {} s".format(time.time() - start_time)
    print "#"*78
    print "Output of 'verdi calculation list':"
    try:
        print subprocess.check_output(
            ["verdi", "calculation", "list"], 
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        print "Note: the command failed, message: {}".format(e.message)

    if calc.has_finished():
        print "Calculation terminated its execution"
        exited_with_timeout = False
        break

if exited_with_timeout:
    print "Timeout!! Calculation did not complete after {} seconds".format(
        timeout_secs)
    sys.exit(2)
else:
    if abs(calc.res.FreeE - expected_free_energy) < 1.e-3:
        print "OK, energy has the expected value"
        sys.exit(0)
    else:
        print "ERROR!"
        print "Expected free energy value: {}".format(expected_free_energy)
        print "Actual free energy value: {}".format(calc.res.FreeE)
        sys.exit(3)
        
